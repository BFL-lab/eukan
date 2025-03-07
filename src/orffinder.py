import re
from Bio.Data import CodonTable
from Bio.Seq import Seq
import pandas as pd
import gffutils
import src.gffparser as gffparser
import src.gffintersecter as gffintersecter


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


def frame_finder(coord):
    if coord % 3 == 1 or coord % 3 == 2:
        return(int(coord % 3 - 1))
    else:
        return(2)


def next_feature_phase_finder(start, end, strand, current_phase):
    if strand == "+":
        return((3 - ((end - start + 1 - current_phase) % 3)) % 3)
    else:
        return((3 - ((start - end + 1 - current_phase) % 3)) % 3)


def consolidate_mRNA_of_overlapping_genes(gff3db, dup_genes):
    """
    Take a gff3db and modify mRNA Parent feature to new gene parent and remove
    the old parent feature fully spanned by the new gene feature
    """
    for f in gff3db.all_features():
        if f.featuretype == "gene" and f.id not in dup_genes:
            yield f
        elif f.featuretype == "mRNA" and f.attributes['Parent'][0] in dup_genes:
            f.attributes['Parent'] = [dup_genes[f.attributes['Parent'][0]]]
            yield f
        elif f.featuretype == "mRNA" and f.attributes['Parent'][0] not in dup_genes:
            yield f
        elif f.featuretype in ['CDS', 'exon']:
            yield f


def fetch_align_fasta(gff3db, fasta, featuretype='exon'):
    seqs = list()
    for mRNA in gff3db.features_of_type('mRNA'):
        seq = str()
        for exon in gff3db.children(mRNA, featuretype=featuretype):
            seq += exon.sequence(fasta, use_strand=False).upper()
        seq = str(Seq(seq).reverse_complement()) if mRNA.strand == '-' else seq
        seqs.append([mRNA.id, mRNA.strand, seq])
    return(seqs)


def find_starts_stops(fastalist, ncbi_genetic_code=1):
    code = CodonTable.unambiguous_dna_by_id[ncbi_genetic_code]
    returnlist = list()
    for name, strand, sequence in fastalist:
        init_term_codons = [
            '(% s)' % '|'.join(code.start_codons),
            '(% s)' % '|'.join(code.stop_codons)]
        seqlen = len(sequence)
        for motif in init_term_codons:
            codon = 'start' if init_term_codons.index(motif) == 0 \
                else 'stop'
            for match in re.finditer(motif, sequence):
                pos = match.start() + 1
                frame = frame_finder(pos)
                returnlist.append(
                    [name, seqlen, pos, strand, codon, frame])
    returnlist = pd.DataFrame(
        returnlist,
        columns=('id', 'seqlen', 'pos', 'strand', 'codon', 'frame'))
    return(returnlist)


def fetch_longest_orf(orf_coords_df, min_orf_len=90, orf_len_frac=0.50):
    # generate combinations of in-phase start+stop combinations of transcript
    b = orf_coords_df[orf_coords_df['codon'] == 'start'].\
        rename(columns={'pos': 'start'}).drop(columns='codon')
    e = orf_coords_df[orf_coords_df['codon'] == 'stop'].\
        rename(columns={'pos': 'stop'}).drop(columns='codon')
    be = pd.merge(b, e, on=['id', 'seqlen', 'strand', 'frame'])
    be = be[be['start'] < be['stop']].reset_index()
    minstop = be.groupby(
            ['id', 'seqlen', 'start', 'strand', 'frame'], as_index=False).\
        apply(lambda x: min(x.stop)).reset_index().\
        rename(columns={None: 'minstop'})
    be = pd.merge(be, minstop,
        left_on = ['id', 'seqlen', 'start', 'strand', 'frame'],
        right_on = ['id', 'seqlen', 'start', 'strand', 'frame'])
    be['orflen'] = be['minstop'] - be['start']
    be = be[be['orflen'] == be['stop'] - be['start']].\
        drop(['minstop', 'index_x'], axis=1)
    # filter for the longest orf
    be = be.loc[be.groupby(['id', 'seqlen', 'strand'])['orflen'].idxmax()]
    # enforce min orf length and orf length fraction
    be = be[
        (be['orflen'] > min_orf_len) &
        (be['orflen']/be['seqlen'] > orf_len_frac)]
    be['stop'] += 2
    return(be)


def orf_to_genome_coords(longest_orfs, gff3):
    for mRNA in gff3.features_of_type('mRNA'):
        if mRNA.id in longest_orfs['id'].values:
            # extract entry from longest_orf
            orf = longest_orfs.loc[longest_orfs.id == mRNA.id]
            # fetch contiguous orf coordinates on transcript, adjust to strand
            orf_coords = (orf.start.item(), orf.stop.item())
            num_exons = len([f.id for f in gff3.children('asmbl_1_mRNA', featuretype='exon')])
            index_shift = 0
            exon_index = 0
            cds_index = 0
            next_feature_phase = 0
            if mRNA.strand == "+":
                exons = gff3.children(mRNA, featuretype="exon", order_by="start")
                for exon in exons:
                    exon_index += 1
                    e_end = exon.end - exon.start
                    prev_orf_coords = orf_coords
                    orf_coords = (
                        max(prev_orf_coords[0] - e_end, 0),
                        orf_coords[1] - e_end - 1)
                    if prev_orf_coords[0] == e_end and e_end > 0:
                        cds_index += 1
                        cds_start = exon.end - exon.start - 1
                        cds_end = exon.end - exon.start
                        index_shift = exon_index
                    elif prev_orf_coords[0] >= e_end:
                        orf_coords = (orf_coords[0], orf_coords[1] + 1)
                        continue
                    elif prev_orf_coords[0] >= 0 \
                            and prev_orf_coords[0] < e_end \
                            and prev_orf_coords[1] > 0:
                        index_shift = exon_index if cds_index == 0 else index_shift
                        cds_index += 1
                        cds_start = max(0, prev_orf_coords[0] - index_shift)
                        if e_end == mRNA.end:
                            cds_end = min(e_end - 1, prev_orf_coords[1])
                        else:
                            cds_end = min(e_end, prev_orf_coords[1] - index_shift)
                        if cds_start > cds_end:
                            continue
                    else:
                        continue
                    if cds_index == 1:
                        phase = 0
                        next_feature_phase = (3 - ((cds_end - cds_start + 1 - phase) % 3)) % 3
                    else:
                        phase = next_feature_phase
                        next_feature_phase = (3 - ((cds_end - cds_start + 1 - phase) % 3)) % 3
                    # create a new CDS feature with coordinates inferred above
                    ID = re.sub('exon', 'CDS', exon.attributes['ID'][0])
                    exon.attributes['ID'] = ID
                    CDS = gffutils.Feature(
                        seqid=exon.chrom, source=exon.source, featuretype='CDS',
                        start=cds_start + exon.start,
                        end=cds_end + exon.start, strand=exon.strand, frame=phase,
                        attributes=exon.attributes)
                    yield CDS
            else:
                exons = gff3.children(mRNA, featuretype="exon", order_by="start", reverse=True)
                orf_coords = (orf_coords[0] - 1, orf_coords[1] - 1)
                for exon in exons:
                    exon_index += 1
                    e_end = exon.end - exon.start
                    prev_orf_coords = orf_coords
                    orf_coords = (
                        max(prev_orf_coords[0] - e_end, 0),
                        orf_coords[1] - e_end - 1)
                    if prev_orf_coords[0] >= e_end:
                        orf_coords = (orf_coords[0] - 1, orf_coords[1])
                        continue
                    elif prev_orf_coords[0] >= 0 \
                            and prev_orf_coords[0] < e_end \
                            and prev_orf_coords[1] > 0:
                        cds_index += 1
                        cds_start = max(0, prev_orf_coords[0])
                        cds_end = min(e_end, prev_orf_coords[1])
                    else:
                        continue
                    if cds_index == 1:
                        phase = 0
                    else:
                        phase = next_feature_phase
                    # create a new CDS feature with coordinates inferred above
                    next_feature_phase = (3 - ((cds_end - cds_start + 1 - phase) % 3)) % 3
                    ID = re.sub('exon', 'CDS', exon.attributes['ID'][0])
                    exon.attributes['ID'] = ID
                    CDS = gffutils.Feature(
                        seqid=exon.chrom, source=exon.source, featuretype='CDS',
                        start=exon.end - cds_end, strand=exon.strand, frame=phase,
                        end=exon.end - cds_start,
                        attributes=exon.attributes)
                    yield CDS


def create_transcriptome_orf_db(gff, genome):
    dialect = gffutils.DataIterator(gff).dialect
    dialect['fmt'] = 'gtf'
    xcripts = gffutils.create_db(
        gff, ':memory:',
        dialect=dialect, gtf_transcript_key='Parent', gtf_gene_key='Parent')
    dialect['fmt'] = 'gff3'
    xcripts.dialect = dialect
    xcripts = gffutils.create_db(
        xcripts, ':memory:', transform=gffparser.gff3_it)
    xcripts_seqs = fetch_align_fasta(xcripts, genome)
    starts_stops = find_starts_stops(xcripts_seqs)
    longest_orfs = fetch_longest_orf(starts_stops)
    # longest orf seems to be okay, maybe translating to genome coords is
    # problem?
    xcripts.update(orf_to_genome_coords(longest_orfs, xcripts))
    consolidated_xcripts = gffintersecter.\
        merge_fully_overlapping_transcript_genes(xcripts)
    return(consolidated_xcripts)
