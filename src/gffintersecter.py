from gffutils.pybedtools_integration import to_bedtool
import pybedtools
import os
import re
import tempfile
import src.gffparser as gffparser
import src.orffinder as orffinder
import gffutils
import pandas as pd


def FeatureDB_to_bedtool_iterator(FeatureDB, features='all'):
    if type(features) is list:
        for feat in features:
            for entry in FeatureDB.features_of_type(feat):
                yield entry
    elif features == 'all':
        for entry in FeatureDB.all_features():
            yield entry
    elif type(features) is str:
        for entry in FeatureDB.features_of_type(features):
            yield entry
    else:
        pass


def merge_fully_overlapping_transcript_genes(gff3db):
    """
    1. find overlapping genes on same strand
    2. take longest of overlapping clusters, associate shorter ones with longest
    3. run through mRNA features, modify parent where there's a match
    """
    pybedtools.set_tempdir(os.getcwd())
    s = to_bedtool(FeatureDB_to_bedtool_iterator(gff3db)).\
        saveas(tempfile.NamedTemporaryFile(dir=os.getcwd()).name)
    genes = s.filter(lambda x: x[2] == "gene").sort()
    ss = genes.merge(s=True, d=-1, c=9, o='distinct_only').\
        to_dataframe(disable_auto_names=True, header=None, low_memory=False)
    list_of_overlaps = [f.replace(";,", ";").replace("ID=", "").split(";")[:-1] for f in ss[3]]
    replacement_ids = {}
    for sublist in list_of_overlaps:
        if len(sublist) > 1:
            for spanning_gene in sublist[1:]:
                replacement_ids[spanning_gene] = sublist[0]
    gff3db_mod = gffutils.create_db(
        orffinder.consolidate_mRNA_of_overlapping_genes(gff3db, replacement_ids),
        ':memory:')
    return(gff3db_mod)


def find_nonoverlapping_genes(gff3_db1, gff3_db2):
    """
    Return all ORF-containing, non-overlapping genes (and their sub-features) from gff3_1 in gff3_2
    """
    pybedtools.set_tempdir(os.getcwd())
    s = to_bedtool(FeatureDB_to_bedtool_iterator(gff3_db1)).\
        saveas(tempfile.NamedTemporaryFile(dir=os.getcwd()).name)
    t = to_bedtool(FeatureDB_to_bedtool_iterator(gff3_db2)).\
        saveas(tempfile.NamedTemporaryFile(dir=os.getcwd()).name)
    st = s.intersect(t, loj=True, s=True).\
        to_dataframe(disable_auto_names=True, header=None, low_memory=False)
    non_overlapping_transcripts = st[(st[2] == "mRNA") & (st[9] == ".")][8]\
        .str.extract(r'(Parent=[A-Za-z_0-9]+;)').\
        replace(r'Parent=', '', regex=True).\
        replace(r';', '', regex = True)[0]
    genes = list(non_overlapping_transcripts)
    genes_with_orfs = list()
    for gene in genes:
        featuretypes = [f.featuretype for f in gff3_db1.children(gene)]
        if 'CDS' in featuretypes:
            genes_with_orfs.append(gene)
    orf_gene_entries = []
    for gene in genes_with_orfs:
        entry = gff3_db1[gene]
        orf_gene_entries.append(entry)
        gene_children = gff3_db1.children(gene, order_by = 'featuretype', reverse = True)
        orf_gene_entries.extend(list(gene_children))
    return(orf_gene_entries)


def find_concordant_models(gff3_1, gff3_2):
    pybedtools.set_tempdir(os.getcwd())
    db1 = gffutils.create_db(gff3_1, ':memory:')
    db2 = gffutils.create_db(gff3_2, ':memory:')
    s = to_bedtool(FeatureDB_to_bedtool_iterator(db1)).\
        saveas(tempfile.NamedTemporaryFile(dir=os.getcwd()).name)
    t = to_bedtool(FeatureDB_to_bedtool_iterator(db2)).\
        saveas(tempfile.NamedTemporaryFile(dir=os.getcwd()).name)
    st = s.intersect(t, wo=True).\
        to_dataframe(disable_auto_names=True, header=None, low_memory=False)
    st = st[st[2] == st[11]]
    # find number of CDS entries per mRNA for both a and b to use as filtering
    sdf = s.to_dataframe(low_memory=False)
    sdf_attrs = attrs_extract(sdf['attributes'])
    sdf = pd.concat([sdf, sdf_attrs], axis=1, sort=False).\
        drop('attributes', axis=1).\
        rename(columns={"ID": "sID", "Parent": "sParent", "Name": "sName"})
    s_numCDS = sdf[sdf['feature'] == 'CDS'].\
        groupby('sParent')['feature'].\
        agg(s_numCDSfeats='count').\
        reset_index()
    tdf = t.to_dataframe(low_memory=False)
    tdf_attrs = attrs_extract(tdf['attributes'])
    tdf = pd.concat([tdf, tdf_attrs], axis=1, sort=False).\
        drop('attributes', axis=1).\
        rename(columns={"ID": "tID", "Parent": "tParent", "Name": "tName"})
    t_numCDS = tdf[tdf['feature'] == 'CDS'].\
        groupby('tParent')['feature'].\
        agg(t_numCDSfeats='count').\
        reset_index()
    # convert attributes to columns
    st_sattrs = attrs_extract(st[8]).add_prefix('s')
    st_tattrs = attrs_extract(st[17]).add_prefix('t')
    stexplicit = pd.concat([st, st_sattrs, st_tattrs], axis=1)
    stexplicit_mRNA = stexplicit[stexplicit[2] == "mRNA"]
    stexplicit_mRNA = pd.merge(
            stexplicit_mRNA, s_numCDS,
            left_on='sID', right_on='sParent', how='left').\
        rename(columns={"sParent_x": "sParent"}).\
        drop(['sParent_y', 8, 17], axis=1)
    stexplicit_mRNA = pd.merge(
        stexplicit_mRNA, t_numCDS,
        left_on='tID', right_on='tParent', how='left').\
        drop('tParent_y', axis=1).\
        rename(columns={'tParent_x': 'tParent'})
    # compute mRNA subject/query lengths
    stexplicit_mRNA['slen'] = abs(stexplicit_mRNA[4]-stexplicit_mRNA[3])
    stexplicit_mRNA['tlen'] = abs(stexplicit_mRNA[13]-stexplicit_mRNA[12])
    stexplicit_CDS = stexplicit[stexplicit[2] == "CDS"]
    stexplicit_CDS = pd.merge(stexplicit_CDS, s_numCDS, how='left').\
        drop([8, 17], axis=1)
    stexplicit_CDS = pd.merge(stexplicit_CDS, t_numCDS, how='left')
    # we're looking for overlaps that share exact intron-CDS structure with
    # some wiggle room at the 5'/3' of the mRNA
    # first filter out mRNA overlaps that don't share the same number of CDS
    numCDS_match = stexplicit_mRNA[
        stexplicit_mRNA['s_numCDSfeats'] == stexplicit_mRNA['t_numCDSfeats']]
    # filter out mRNA subject overlaps that exceed 5% of the subject
    # length at either the 5' or the 3'
    numCDS_match = numCDS_match[
        (numCDS_match[12] <= numCDS_match[3] + numCDS_match[3]*0.05) &
        (numCDS_match[12] >= numCDS_match[3] - numCDS_match[3]*0.05) &
        (numCDS_match[13] <= numCDS_match[4] + numCDS_match[4]*0.05) &
        (numCDS_match[13] >= numCDS_match[4] - numCDS_match[4]*0.05)]
    # filter out all CDS whose parent-pairs are not in the resulting
    # list from the previous two steps
    matching_mRNA = numCDS_match[['sID', 'tID']].\
        rename(columns={'sID': 'sParent', 'tID': 'tParent'})
    matching_CDS = pd.merge(stexplicit_CDS, matching_mRNA)
    # filter out subfeatures whose parents don't have the same number of
    # pairwise overlapping subfeatures
    matching_CDS['rank'] = matching_CDS.groupby(['sParent', 'tParent']).\
        cumcount()+1
    num_overlapping_feats = matching_CDS.groupby(['sParent', 'tParent']).\
        size().reset_index(name='numfeatoverlap')
    matching_CDS = pd.merge(matching_CDS, num_overlapping_feats, how='left')
    matching_CDS = matching_CDS[
        matching_CDS['numfeatoverlap'] == matching_CDS['s_numCDSfeats']]
    # compute intron junctions (where applicable) and assert junction matching
    # logic
    matching_CDS = matching_CDS.\
        apply(lambda x: match_exon_intron_bounds(x), axis=1)
    matching_CDS = matching_CDS.groupby(['sParent', 'tParent']).\
        filter(
            lambda x:
            len(set(x['bounds'])) == 1 and 'mismatch' not in x['bounds'])
    # now fetch all matching parent feature pairs to recombine the filtered
    # models
    corr_mRNA_pairs = matching_CDS[['sParent', 'tParent']].\
        rename(columns={'sParent': 'sID', 'tParent': 'tID'})
    matching_mRNA = pd.merge(stexplicit_mRNA, corr_mRNA_pairs, how='inner')
    corr_gene_pairs = matching_mRNA[['sParent', 'tParent']].\
        drop_duplicates().\
        rename(columns={'sParent': 'sID', 'tParent': 'tID'})
    matching_gene = pd.merge(stexplicit, corr_gene_pairs, how='inner')
    matching_gene[8] = "ID=" + matching_gene['sID'] + ';'
    matching_gene = matching_gene[[0, 1, 2, 3, 4, 5, 6, 7, 8]].\
        drop_duplicates()
    matching_mRNA[8] = "ID=" + matching_mRNA['sID'] + \
        ";Parent=" + matching_mRNA['sParent'] + ";"
    matching_mRNA = matching_mRNA[[0, 1, 2, 3, 4, 5, 6, 7, 8]].\
        drop_duplicates()
    corr_mRNA_pairs.rename(
        columns={'sID': 'sParent', 'tID': 'tParent'}, inplace=True)
    matching_CDS = pd.merge(stexplicit_CDS, corr_mRNA_pairs, how='inner')
    matching_CDS[8] = "ID=" + matching_CDS['sID'] + ";Parent=" + \
        matching_CDS['sParent'] + ";"
    matching_CDS = matching_CDS[[0, 1, 2, 3, 4, 5, 6, 7, 8]].drop_duplicates()
    exons = stexplicit[(stexplicit[2] == 'exon') & (stexplicit[11] == 'exon')]
    matching_exons = pd.merge(exons, corr_mRNA_pairs)
    matching_exons = matching_exons[[0, 1, 2, 3, 4, 5, 6, 7, 8]].\
        drop_duplicates()
    concordantmodels = pd.concat(
        [matching_gene, matching_mRNA, matching_CDS, matching_exons], axis=0)
    concordantmodels = gffutils.create_db(
        concordantmodels.
        to_csv(None, sep='\t', header=False, index=False),
        ':memory:', from_string=True)
    return(concordantmodels)


def combine_nonredundant_models(*FeatureDBs):
    dbs = [
        to_bedtool(FeatureDB_to_bedtool_iterator(f, 'gene')).
        saveas(tempfile.NamedTemporaryFile(dir=os.getcwd()).name) for f in FeatureDBs]
    if len(dbs) == 2:
        nr = [re.sub('(ID=)|;', '', f[8]) for f in dbs[0] + dbs[1]]
        nr = gffutils.create_db(
            gffparser.filter_FeatureDB_by_model(FeatureDBs[0], nr), ':memory:')
    elif len(dbs) == 3:
        nr = [re.sub('(ID=)|;', '', f[8]) for f in dbs[0] + dbs[1] + dbs[2]]
        p1 = [re.sub('(ID=)|;', '', f[8]) for f in dbs[0] + dbs[1] - dbs[2]]
        p2 = [re.sub('(ID=)|;', '', f[8]) for f in dbs[0] - dbs[1] + dbs[2]]
        p3 = [re.sub('(ID=)|;', '', f[8]) for f in dbs[2] + dbs[1] - dbs[0]]
        nr += p1 + p2
        nr = gffutils.create_db(
            gffparser.filter_FeatureDB_by_model(FeatureDBs[0], nr), ':memory:')
        if p3:
            additional = gffutils.create_db(
                gffparser.filter_FeatureDB_by_model(
                    FeatureDBs[2], p3), ':memory:')
            nr.update(additional)
    return(nr)


def extract_supported_models(*gff3):
    # find single best gene feature match between genemark and proteins
    if len(gff3[0]) == 2:
        training_set = find_concordant_models(gff3[0][0], gff3[0][1])
    # find matching alignments between transcriptome and proteins
    elif len(gff3[0]) == 3:
        pair1 = find_concordant_models(gff3[0][0], gff3[0][1])
        pair2 = find_concordant_models(gff3[0][0], gff3[0][2])
        pair3 = find_concordant_models(gff3[0][2], gff3[0][1])
        training_set = combine_nonredundant_models(pair1, pair2, pair3)
    gffparser.featuredb2gff3_file(training_set, 'training_set.gff3')
    training_set = gffutils.create_db('training_set.gff3', ':memory:')
    return(training_set)


def attrs_extract(attributes_column):
    attrs = attributes_column.str.extractall('([^;]+?)=(?=([^;]*))?').\
        reset_index(level=1, drop=True).\
        set_index(0, append=True)[1].\
        unstack(level=1)
    return(attrs)


def match_exon_intron_bounds(line):
    if line['rank'] == 1 and line['s_numCDSfeats'] != 1:
        # first feature of >=2 in order of ascending coordinates
        if line[13] <= line[4] + 3 and line[13] >= line[4] - 3:
            line['bounds'] = 'match'
            return(line)
        else:
            line['bounds'] = 'mismatch'
            return(line)
    elif line['rank'] == line['s_numCDSfeats'] and line['s_numCDSfeats'] != 1:
        # last feature of >=2 in order of ascending coordinates
        if line[12] <= line[3] + 3 and line[12] >= line[3] - 3:
            line['bounds'] = 'match'
            return(line)
        else:
            line['bounds'] = 'mismatch'
            return(line)
    elif line['s_numCDSfeats'] == 1:
        # single exon feature that already "passes" as matching model
        line['bounds'] = 'match'
        return(line)
    elif line[12] <= line[3] + 3 and \
            line[12] >= line[3] - 3 and \
            line[13] <= line[4] + 3 and \
            line[13] >= line[4] - 3:
        # internal feature of >=2 exons gene model
        line['bounds'] = 'match'
        return(line)
    else:
        line['bounds'] = 'mismatch'
        return(line)
