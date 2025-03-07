import gffutils
import re
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict


def gtf2gff3(f):
    """
    Convert a gtf gffutils.FeatureDB to gff3
    """
    feat = f.featuretype
    if feat in ['gene', 'mRNA', 'transcript', 'exon', 'CDS']:
        if feat == 'gene':
            f.attributes['ID'] = f.attributes.pop('gene_id', None)
        elif feat == 'transcript' or feat == 'mRNA':
            f.attributes['ID'] = f.attributes.pop('transcript_id', None)
            f.attributes['Parent'] = f.attributes.pop('gene_id', None)
            f.featuretype = 'mRNA'
        elif feat == 'exon' or feat == 'CDS':
            f.attributes['Parent'] = f.attributes.pop('transcript_id', None)
            f.attributes.pop('gene_id', None)
            f.attributes['ID'] = f.id
        return(f)


def filter_non_gtf_features(f):
    """
    Remove additional feature lines that genemark outputs as of v4.6+
    """
    if f.featuretype in ['CDS', 'exon']:
        print(f.attributes)
        return(f)


def gff3_to_fasta(gff3, genome, extract_to = 'protein', ncbi_gen_code = 1):
    """
    Fetch CDS coordinates from mRNA features and then extract the corresponding
    sequence from the assembly

    extract_to is either protein or transcript
    ncbi_gen_code is integer for genetic code table
    """
    gff3db = gffutils.create_db(gff3, ':memory:')
    contigs = SeqIO.to_dict(SeqIO.parse(genome, 'fasta'))
    extracted_seqs = []
    for f in gff3db.features_of_type('mRNA'):
        seq = []
        children = gff3db.children(f, featuretype = 'CDS', order_by='start')
        for child in children:
            seq.append(contigs[child.chrom][child.start-1:child.end].seq)
        if seq:
            seq_concat = Seq("".join([str(cds) for cds in seq]))
            if f.strand == "-":
                seq_concat = seq_concat.reverse_complement()
            extracted_seqs.append(SeqRecord(seq_concat, description = "", id = f.id))
    if extract_to == 'protein':
        for seq in extracted_seqs:
            trans = seq.translate(table = CodonTable.unambiguous_dna_by_id[ncbi_gen_code])
            trans.id = f.id
            trans.description = ""
            print(trans.format('fasta'))
    else:
        for seq in extracted_seqs:
            print(seq.format('fasta'))


def add_missing_feats_to_gff3(gff3):
    """
    output from various gene prediction/alignment software produce a variety
    of gff3 flavours, often times incomplete. we need to fix them.
    """
    feats = [f for f in gff3.featuretypes()]
    # complete gene/CDS gffs
    if 'mRNA' not in feats and 'gene' in feats:
        for gene in gff3.features_of_type('gene'):
            attrs = {}
            attrs['ID'] = '%s_mRNA' % gene.attributes['ID'][0]
            attrs['Parent'] = gene.attributes['ID'][0]
            mRNA = gffutils.Feature(
                seqid=gene.chrom,
                source=gene.source,
                featuretype='mRNA',
                start=gene.start,
                end=gene.end,
                strand=gene.strand,
                frame='.',
                attributes=attrs)
            yield mRNA
            gene.attributes['ID'][0] += "_gene"
    # complete gene/mRNA/exon gffs
    if 'CDS' not in feats and 'exon' in feats:
        for exon in gff3.features_of_type('exon'):
            phase = (exon.start-1) % 3 if exon.strand == '-' \
                else (exon.end+1) % 3
            ID = re.sub('exon', 'CDS', exon.id)
            exon.attributes['ID'] = ID
            # create a new CDS feature from exons
            CDS = gffutils.Feature(
                seqid=exon.chrom,
                source=exon.source,
                featuretype='CDS',
                start=exon.start,
                end=exon.end,
                strand=exon.strand,
                frame=phase,
                attributes=exon.attributes)
            yield CDS
    if 'exon' not in feats and 'CDS' in feats:
        for CDS in gff3.features_of_type('CDS'):
            # create a new CDS feature from CDSs
            ID = re.sub('CDS|cds', 'exon', CDS.attributes['ID'][0])
            CDS.attributes['ID'] = ID
            exon = gffutils.Feature(
                seqid=CDS.chrom,
                source=CDS.source,
                featuretype='exon',
                start=CDS.start,
                end=CDS.end,
                strand=CDS.strand,
                frame='.',
                attributes=CDS.attributes)
            yield exon


def add_cq_mRNA(f):
    """
    specifically to fix the incomplete f output of CQ
    """
    f.source = 'codingquarry'
    if f.featuretype == 'CDS':
        f.attributes['Parent'][0] += "_mRNA"
    return f


def fix_dup_IDs(f):
    if f.featuretype == 'mRNA':
        f.id = f.attributes['ID'][0]
    else:
        f.attributes['ID'] = f.id
    return(f)


def fix_spaln_cds_featuretype(f):
    if f.featuretype == 'cds':
        f.attributes['ID'][0] = f.attributes['Parent'][0] + ':' + f.attributes['ID'][0]
        f.featuretype = 'CDS'
    return(f)


def fix_spaln_ids(f):
    if f.featuretype == 'CDS' or f.featuretype == 'exon':
        f.attributes['ID'][0] = f.id
    f.source = 'prot_align'
    return(f)


def prot2augustus_hints(f):
    if f.featuretype == 'CDS':
        group = f.attributes['Parent']
        f.attributes = {'pri': '1', 'src': 'P', 'group': group}
        f.featuretype = 'CDSpart'
        return(f)


def fix_contig_names(f):
    f.chrom = f.chrom.split(' ')[0]
    return(f)


def featuredb2gff3_file(featuredb, out):
    with open(out, 'w') as fout:
        for gene in featuredb.features_of_type('gene', order_by=('seqid', 'start')):
            fout.write('%s\n' % str(gene))
            for mRNA in featuredb.children(
                    gene,
                    featuretype='mRNA',
                    order_by='start'):
                fout.write('%s\n' % str(mRNA))
                for gchild in ['exon', 'CDS']:
                    for f in featuredb.children(
                            mRNA,
                            featuretype=gchild,
                            order_by='start'):
                        fout.write('%s\n' % str(f))


def gff3_it(f):
    if f.featuretype == 'transcript':
        f.featuretype = 'mRNA'
        f.attributes['ID'] = f.attributes['Parent'][0] + '_mRNA'
    elif f.featuretype == 'gene':
        f.attributes['ID'] = f.attributes.pop('Parent', None)
    else:
        f.attributes['Parent'] = f.attributes['Parent'][0] + '_mRNA'
    if not f.attributes['ID'] and f.id:
        f.attributes['ID'] = f.id
    return(f)


def filter_FeatureDB_by_model(featuredb, models):
    for model in models:
        yield featuredb[model]
        for child in featuredb.children(model):
            yield child


def gff2zff(gff3, fasta):
    featuredb = gffutils.create_db(gff3, ':memory:')
    contig_list = [f.id for f in SeqIO.parse(fasta, 'fasta')]
    contigs = defaultdict(list, { k:[] for k in contig_list })
    for mRNA in featuredb.features_of_type('mRNA'):
        tot_exons = len([f.id for f in featuredb.children(mRNA, featuretype='exon')])
        exon_num = 1
        for exon in featuredb.children(mRNA, featuretype='exon'):
            if tot_exons == 1:
                feature = 'Esngl'
            elif exon_num == 1 and tot_exons > 1:
                feature = 'Einit' if mRNA.strand == '+' else 'Eterm'
            elif exon_num > 1 and exon_num < tot_exons:
                feature = 'Exon'
            elif exon_num == tot_exons:
                feature = 'Eterm' if mRNA.strand == '+' else 'Einit'
            coords = (exon.start, exon.end) \
                if mRNA.strand == '+' \
                else (exon.end, exon.start)
            parent = exon.attributes['Parent'][0]
            contigs[mRNA.chrom].append([feature, coords[0], coords[1], parent])
            exon_num += 1
    with open('genome.ann', 'w') as fout:
        for contig in contigs:
            fout.write('>%s\n' % contig)
            for i in contigs[contig]:
                fout.write("%s  %s %s %s\n" % (i[0], i[1], i[2], i[3]))


def fix_snap_featuretype(f):
    f.source = 'snap'
    f.featuretype = 'exon'
    f.attributes['Parent'] = list(f.attributes.keys())[0]
    f.attributes['ID'] = []
    f.attributes = {'ID': f.attributes['ID'], 'Parent': f.attributes['Parent']}
    return(f)


def homogenize_snap_source(f):
    f.source = 'snap'
    return(f)


def protgff3_to_evm(gff3):
    if gff3.featuretype == 'CDS':
        gff3.source = 'spaln'
        gff3.featuretype = 'nucleotide_to_protein_match'
        return(gff3)


def clean_augustus_gff3(f):
    if f.featuretype in ['gene', 'mRNA', 'CDS', 'exon']:
        f.source = 'augustus'
        return(f)


def filter_gth_gff3(f):
    if f.featuretype in ['gene', 'exon']:
        if f.featuretype == 'exon':
            f.attributes['Parent'][0] += '_mRNA'
        return(f)


def fix_gth_ids(f):
    f.source = 'prot_align'
    if f.featuretype == 'exon':
        f.attributes['ID'] = [f.id]
    if f.featuretype == 'mRNA':
        f.id = f.attributes['ID'][0]
    return(f)


def prettify_gff3(gff3, shortname):
    gene_num = 1
    for gene in gff3.features_of_type('gene', order_by=('seqid', 'start')):
        gene_name = '%s_%s' % (shortname, str(gene_num).zfill(5))
        gene.attributes['ID'] = gene_name
        gene.attributes['locus_tag'] = gene_name
        gene.attributes['Name'] = gene_name
        gene.source = 'eukannotpass'
        gene_num += 1
        mRNA_num = 1
        yield gene
        for mRNA in gff3.children(gene, featuretype='mRNA', order_by=('seqid', 'start')):
            mRNA_name = gene_name + '.mRNA.' + str(mRNA_num)
            mRNA.attributes['ID'] = mRNA_name
            mRNA.attributes['Name'] = mRNA_name
            mRNA.attributes['Parent'] = gene_name
            mRNA.attributes['locus_tag'] = mRNA_name
            mRNA.source = 'eukannotpass'
            mRNA_num += 1
            yield mRNA
            for gchild in list(set([f.featuretype for f in gff3.children(gene, level=2)])):
                f_num = 1
                for f in gff3.children(mRNA, featuretype=gchild):
                    f_name = '%s:%s:%s' % (mRNA_name, f.featuretype, str(f_num))
                    f.attributes['ID'] = f_name
                    f.attributes['Parent'] = mRNA_name
                    f.attributes['locus_tag'] = f_name
                    f.source = 'eukannotpass'
                    f_num += 1
                    yield f


def fix_CDS_phases(gff3):
    for gene in gff3.features_of_type('gene'):
        yield gene
        for mRNA in gff3.children(gene, featuretype="mRNA"):
            cds_index = 0
            next_feature_phase = 0
            yield mRNA
            for exon in gff3.children(mRNA, featuretype='exon'):
                yield exon
            if mRNA.strand == "+":
                num_cds = len([f for f in gff3.children(mRNA, featuretype='CDS')])
                if num_cds == 0:
                    continue
                all_cds = gff3.children(mRNA, featuretype="CDS", order_by="start")
                for cds in all_cds:
                    cds_index += 1
                    if cds_index == 1:
                        phase = 0
                        next_feature_phase = (3 - ((cds.end - cds.start + 1 - phase) % 3)) % 3
                    else:
                        phase = next_feature_phase
                        next_feature_phase = (3 - ((cds.end - cds.start + 1 - phase) % 3)) % 3
                    # create a new CDS feature with coordinates inferred above
                    cds.frame = str(phase)
                    yield cds
            else:
                all_cds = gff3.children(mRNA, featuretype="CDS", order_by="start", reverse=True)
                for cds in all_cds:
                    cds_index += 1
                    if cds_index == 1:
                        phase = 0
                    else:
                        phase = next_feature_phase
                    # create a new CDS feature with coordinates inferred above
                    next_feature_phase = (3 - ((cds.end - cds.start + 1 - phase) % 3)) % 3
                    cds.frame = str(phase)
                    yield cds
