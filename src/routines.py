import os
import re
import sys
import glob
import shutil
import string
import datetime
import gffutils
import src.gffparser as gffparser
import src.gffintersecter as gffintersecter
import src.orffinder as orffinder
import multiprocessing as mp
import pandas as pd
import random
import subprocess
import statistics
from Bio import SeqIO
from Bio.Data import CodonTable
from fileinput import FileInput
from pathlib import Path


class Routines:
    """
    Holds the routines for repeat library construction, spliced alignment of
    proteins, gene predictions and consensus model creation
    """
    def __init__(self, args):
        args.genome = os.path.abspath(args.genome)
        is_fasta(args.genome)
        self.wd = os.getcwd()
        self.name = os.path.splitext(os.path.basename(args.genome))[0]
        self.shortname = _rand_string(5)
        self.gff3_dialect = {
            'leading semicolon': False,
            'trailing semicolon': True,
            'quoted GFF2 values': False, 'field separator': ';',
            'keyval separator': '=', 'multival separator': ',', 'fmt': 'gff3',
            'repeated keys': False, 'order': ['ID', 'Name', 'Parent']}
        self.args = args
        self.numcpu = int(self.args.numcpu)

    def find_orfs(self, trans_gff3):
        # infer orfs in the transcripts
        is_gff(trans_gff3)
        if _init('orf_finder', 'transcript_orfs.gff3'):
            return('%s/orf_finder/transcript_orfs.gff3' % self.wd)
        orfs = orffinder.create_transcriptome_orf_db(
            trans_gff3, self.args.genome)
        gffparser.featuredb2gff3_file(orfs, 'transcript_orfs.gff3')
        xcript_orfs = '%s/orf_finder/transcript_orfs.gff3' % self.wd
        os.chdir(self.wd)
        sys.stdout.write(' done.\n')
        return(xcript_orfs)

    def align_proteins(self, gff3, proteins, *evidence):
        """
        Spliced alignment of protein sequences, specified on the command line,
        to the genome
        """
        if _init('prot_align', 'prot.gff3'):
            return('%s/prot_align/prot.gff3' % self.wd)
        with open('prots.faa', 'w') as outfile:
            for filename in proteins:
                with open(filename) as infile:
                    outfile.write(infile.read())
        is_fasta('prots.faa')
        os.symlink(self.args.genome, "genome")
        os.symlink(self.args.genome, "genome.gf")
        models = gffutils.create_db(
            gff3, ':memory:', merge_strategy='create_unique')
        gene_count = 0
        intron_count = 0
        for f in models.features_of_type('gene'):
            gene_count += 1
            introns_in_gene = len(
                [e for e in models.children(f, featuretype='exon')]) - 1
            intron_count += introns_in_gene
        if intron_count/gene_count > 0.25:  # run spaln
            lengths = [
                abs(f.end - f.start) for f in models.features_of_type('gene')]
            mean_len = int(sum(lengths)/len(lengths))
            sd_len = int(statistics.stdev(lengths))
            maxlen = mean_len + 2*sd_len
            spaln_parms = ['-XG%d' % maxlen, '-C70', '-C30']
            # build a genome index for protein alignments
            _run("makdbs -KD genome")
            # make the associated index files for alignment
            _run("makblk.pl -Wgenome.bkp genome.gf")
            # compute intron length distribution
            intron_hints = evidence[0] if evidence else gff3
            introns = pd.read_csv(intron_hints, sep='\t', header=None)
            introns = pd.DataFrame((introns[4] - introns[3]).value_counts()) \
                .reset_index() \
                .sort_values(by='index')
            introns.to_csv('introns.ild', sep=' ', header=None, index=False)
            # infer intron length distribution parameters
            int_dist = _run("fitild -d IldModel.txt introns.ild")
            with open("cmd.out", "r") as fp:
                for line in fp:
                    if re.search("^introns\s", line):
                        int_dist = [
                            "{:f}".format(float(x)) for x in
                            line.replace('\t', ' ').
                            replace('\n', ' ').
                            split()[3:-3]]
            del int_dist[1]
            int_dist = ' '.join(str(x) for x in int_dist)
            # modify alignment parameters in the index
            _run("spaln -Wgenome.bkp %s -KP genome.gf" % spaln_parms[0])
            # run the spliced alignment with the above parameters
            cmd = "nice -19 spaln -O12 -Q7 -yX -LS " + \
                "-C%s -t%s -yI\"%s\" -dgenome prots.faa" % \
                (self.args.code, self.numcpu, int_dist)
            with open('cmd.log', 'a') as f:
                f.write("%s\ncmd: %s\n" % (datetime.datetime.now(), cmd))
            p = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            with open('cmd.log', 'a') as f:
                f.write('%s\n' % cmd)
            # convert alignments to gff format
            _run('sortgrcd -F0 -C60 -O0 prots.grd', 'temp.gff3')
            spaln_out = gffutils.create_db(
                'temp.gff3', ':memory:',
                transform=gffparser.fix_spaln_cds_featuretype,
                merge_strategy='create_unique')
            spaln_out = gffutils.create_db(
                spaln_out, ':memory:',
                transform=gffparser.fix_spaln_ids)
            # add missing exon features and correct improperly formatted cds
            spaln_out.update(gffparser.add_missing_feats_to_gff3(spaln_out))
            gffparser.featuredb2gff3_file(spaln_out, 'prot.gff3')
        else:  # run gth
            _run(
                "gth -genomic genome.gf " +
                "-protein prots.faa -intermediate -xmlout -gzip " +
                "-o gth.xml.gz -force")
            _run(
                "gthconsensus -mincoverage 0.70 -minalignmentscore 0.30 " +
                "-intermediate -gff3out -force -o gth.gff3 gth.xml.gz")
            gth_out = gffutils.create_db(
                'gth.gff3', ':memory:', transform=gffparser.filter_gth_gff3,
                merge_strategy='create_unique')
            gth_out.update(
                gffparser.add_missing_feats_to_gff3(gth_out),
                merge_strategy='create_unique')
            gth_out = gffutils.create_db(
                gth_out, ':memory:', transform=gffparser.fix_gth_ids)
            gffparser.featuredb2gff3_file(gth_out, 'prot.gff3')
        # create a CDS hints file from alignments
        hints = gffutils.create_db(
            'prot.gff3', ':memory:',
            transform=gffparser.prot2augustus_hints)
        with open('hints_protein.gff', 'w') as fout:
            for feat in hints.all_features():
                fout.write('%s\n' % str(feat))
        protgff3 = '%s/prot_align/prot.gff3' % self.wd
        os.chdir(self.wd)
        sys.stdout.write(' done.\n')
        return(protgff3)

    def run_genemark(self, *evidence):
        """
        Subroutine to build and run a GeneMark.HMM command
        """
        if self.args.rnaseq_hints is None:
            is_gff(self.args.rnaseq_hints)
        if _init('genemark', 'genemark.gff3'):
            return('%s/genemark/genemark.gff3' % self.wd)
        # turn on high gene density switch if fungus
        hgdswitch = '--fungus' if self.args.fungus is True else ''
        # if no RNA-Seq available, run genemark in self-training mode
        if self.args.rnaseq_hints is None:
            trainingtype = '--ES'
        else:
            # else, extract the intron data for guided training
            hints = pd.read_csv(
                evidence[0], sep="\t", header=None, low_memory=False)
            introns = hints[hints[2] == 'intron']
            introns.to_csv('introns.gff', sep='\t', header=None, index=False)
            # check if number of introns is >=150, otherwise GM breaks on ET
            trainingtype = '--ET=introns.gff --et_score=3' \
                if introns.shape[0] >= 150 else '--ES'
        # run genemark
        _run(
            "nice -19 gmes_petap.pl --soft 1000 " +
            "%s --cores=%s --sequence=%s %s"
            % (trainingtype, self.numcpu, self.args.genome, hgdswitch))
        # load in the genemark gtf
        gmgtf = gffutils.create_db('genemark.gtf', ':memory:', verbose=False)
        gmgff3 = gffutils.create_db(
            gmgtf, ':memory:', dialect=self.gff3_dialect,
            transform=gffparser.gtf2gff3, verbose=False)
        gmgff3.update(
            gffparser.add_missing_feats_to_gff3(gmgff3),
            merge_strategy='create_unique')
        gmgff3 = gffutils.create_db(gmgff3, ':memory:',
            transform=gffparser.fix_contig_names, verbose=False)
        gffparser.featuredb2gff3_file(gmgff3, 'genemark.gff3')
        gmgff3 = '%s/genemark/genemark.gff3' % self.wd
        os.chdir(self.wd)
        sys.stdout.write(' done.\n')
        return(gmgff3)

    def run_codingquarry(self, evidence):
        """
        Subroutine to build and run a CodingQuarry commands
        """
        if _init('codingquarry', 'codingquarry.gff3'):
            return('%s/codingquarry/codingquarry.gff3' % self.wd)
        # use the rna-seq assembly if it's available
        flag = '-t' if self.args.transcriptsGFF else '-a'
        hints = '%s %s' % (flag, evidence)
        # run codingquarry
        _run(
            "nice -19 CodingQuarry -p %s -f %s %s" %
            (self.numcpu, self.args.genome, hints))
        cq = gffutils.create_db(
            'out/PredictedPass.gff3', ':memory:',
            transform=gffparser.add_cq_mRNA,
            merge_strategy='create_unique')
        cq.update(
            gffparser.add_missing_feats_to_gff3(cq),
            merge_strategy='create_unique')
        cq = gffutils.create_db(
            cq, ':memory:', transform=gffparser.fix_dup_IDs,
            merge_strategy='create_unique')
        gffparser.featuredb2gff3_file(cq, 'codingquarry.gff3')
        cqgff3 = '%s/codingquarry/codingquarry.gff3' % self.wd
        os.chdir(self.wd)
        sys.stdout.write(' done.\n')
        return(cqgff3)

    def build_training_set(self, *gff3):
        """
        Subroutine to build evidence-driven training set for gene predictors
        """
        modelsdb = gffintersecter.extract_supported_models(gff3[0])
        # find single best gene feature match between genemark and proteins
        num_training_models = int(len(
            [f.id for f in modelsdb.features_of_type('gene')])*1/4)
        # extract a training set from the genemark models
        gff3db = gffutils.create_db(gff3[0][0], ':memory:')
        gene_lengths = [
            f.end - f.start for f in gff3db.features_of_type('gene')]
        flank = int(sum(gene_lengths)/len(gene_lengths)/2)
        _run(
            "gff2gbSmallDNA.pl training_set.gff3 " +
            "%s %s genbank.gb" % (self.args.genome, flank))
        _run("randomSplit.pl genbank.gb %d" % num_training_models)
        return

    def run_augustus(self, *evidence):
        """
        Subroutine to run augustus commands
        """
        if _init('augustus', 'augustus.gff3'):
            return('%s/augustus/augustus.gff3' % self.wd)
        # gather all hints together
        filenames = ['%s/prot_align/hints_protein.gff' % self.wd]
        if self.args.rnaseq_hints:
            ext_cfg = 'eukka.MPEW.RM.cfg'
            filenames.append('%s/%s' % (self.wd, self.args.rnaseq_hints))
        else:
            ext_cfg = "eukka.MP.cfg"
        # define and create a new species to store augustus parameters
        config_path = os.environ['AUGUSTUS_CONFIG_PATH']
        aug_config = "%s/species/%s" % (config_path, self.name)
        if os.path.exists(aug_config):
            shutil.rmtree(aug_config)
        _run("new_species.pl --species=%s" % self.name)
        # extract a training set from the genemark models
        self.build_training_set(evidence)
        # first round of training
        _run("etraining --species=%s genbank.gb.train" % self.name)
        # modify the stop codon frequency based on training set
        if self.args.code == '6':
            stopprob = ['0', '0', '1']
        else:
            stopprob = list()
            with open("cmd.out", "r") as fp:
                for line in fp:
                    if re.search("(tag:|taa:|tga:)", line):
                        stopprob.append(re.sub('[()]', '', line).split()[-1])
        aug_config_file = "%s/species/%s/%s_parameters.cfg" % \
            (os.environ['AUGUSTUS_CONFIG_PATH'], self.name, self.name)
        with FileInput(files=[aug_config_file], inplace=True) as f:
            for line in f:
                line = line.rstrip()
                if re.search("(amberprob|ochreprob|opalprob)", line):
                    stopline = line.split()
                    stopline[1] = stopprob.pop(0)
                    stopline = ' '.join(stopline)
                    print(stopline)
                else:
                    print(line)
        # run predictor on the test models
        _run("augustus --species=%s genbank.gb.test" % self.name)
        # optimize the parameters
        _run(
            "nice -19 optimize_augustus.pl --species=%s " % self.name +
            "--onlytrain=genbank.gb.train --cpus=%s " % self.numcpu +
            "genbank.gb.test")
        # second round
        _run("etraining --species=%s genbank.gb.train " % self.name)
        _run("augustus --species=%s genbank.gb.test " % self.name)
        # create the hints file
        with open('hints_all.gff', 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
        # make the final predictions; split the genome and run in parallel
        handle = open(self.args.genome, "rU")
        assemblysize = int()
        for record in SeqIO.parse(handle, "fasta"):
            assemblysize += len(record)
        handle.close()
        minsize = int(assemblysize/self.numcpu*2)
        os.symlink(self.args.genome, 'genome.fa')
        _run("splitMfasta.pl genome.fa --minsize=%s" % minsize)
        splits = splits = glob.glob('./genome.split.*.fa')
        base_cmd = 'nice -19 augustus --species=%s ' % self.name + \
            '--extrinsicCfgFile=%s/extrinsic/%s ' % (config_path, ext_cfg) + \
            '--hintsfile=hints_all.gff --softmasking=1 --UTR=off '
        cmds = [
            mp.Process(target=_run, args=(base_cmd + split, '%s.gff3' % split))
            for split in splits]
        cmds_start = [cmd.start() for cmd in cmds]
        cmds_wait = [cmd.join() for cmd in cmds]
        f = open('augustus.gff', 'w')
        join = subprocess.Popen(
            'cat genome.split.*gff3 | join_aug_pred.pl | gtf2gff.pl ' +
            '--printExon -gff3 --out=augustus.gff', shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        join_out = join.communicate()[0].decode('utf-8')
        f.write(join_out)
        f.close()
        aug_out = gffutils.create_db(
            'augustus.gff', ':memory:',
            transform=gffparser.clean_augustus_gff3)
        gffparser.featuredb2gff3_file(aug_out, 'augustus.gff3')
        map(os.remove, glob.glob("*split*"))
        auggff3 = '%s/augustus/augustus.gff3' % self.wd
        os.chdir(self.wd)
        sys.stdout.write(' done.\n')
        return(auggff3)

    def run_snap(self, *evidence):
        """
        Subroutine to train and run SNAP
        """
        if _init('snap', 'snap.gff3'):
            return('%s/snap/snap.gff3' % self.wd)
        os.symlink(self.args.genome, 'genome.dna')
        # pull out the best augustus models for training
        self.build_training_set(evidence)
        # convert the gff3 into snap's silly format
        gffparser.gff2zff('training_set.gff3', self.args.genome)
        # construct the snap HMMs
        _run("fathom -categorize 1000 genome.ann genome.dna")
        _run("fathom -export 1000 -plus uni.ann uni.dna")
        _run("forge export.ann export.dna")
        _run('hmm-assembler.pl snap .',  'snap.hmm')
        # run snap using the HMM
        _run(
            "nice -19 snap -gff snap.hmm %s " % self.args.genome, 'snap.gff')
        snap = gffutils.create_db(
            'snap.gff', ':memory:', transform=gffparser.fix_snap_featuretype)
        dialect = gffutils.DataIterator('snap.gff').dialect
        dialect['fmt'] = 'gtf'
        snap = gffutils.create_db(
            snap, ':memory:', dialect=dialect,
            gtf_transcript_key='Parent', gtf_gene_key='Parent')
        snap.dialect = self.gff3_dialect
        snap = gffutils.create_db(
            snap, ':memory:', transform=gffparser.gff3_it)
        snap.update(gffparser.add_missing_feats_to_gff3(snap))
        snap = gffutils.create_db(
            snap, ':memory:', transform=gffparser.homogenize_snap_source)
        gffparser.featuredb2gff3_file(snap, 'snap.gff3')
        snapgff3 = '%s/snap/snap.gff3' % self.wd
        os.chdir(self.wd)
        sys.stdout.write(' done.\n')
        return(snapgff3)

    def run_evm(self, *evidence):
        """
        Routine to run EVidenceModeler with a predefined set of weights
        """
        if _init('evm_consensus_models', '%s.gff3' % self.name):
            return('%s/evm_consensus_models/%s.gff3' % (self.wd, self.name))
        _run("cdbfasta %s" % self.args.genome)
        w = [str(d) for d in self.args.weights]
        weights = {
            'prot.gff3': ['PROTEIN', 'ALN', w[0]],
            'augustus.gff3': ['ABINITIO_PREDICTION', 'augustus', w[1]],
            'snap.gff3': ['ABINITIO_PREDICTION', 'snap', w[1]],
            'genemark.gff3': ['ABINITIO_PREDICTION', 'genemark', w[1]],
            'codingquarry.gff3': ['ABINITIO_PREDICTION', 'codingquarry', w[1]],
            'nr_transcripts.gff3': ['TRANSCRIPT', 'PASA-assembly', w[2]]}
        evm_opts = str()
        with open("weights.txt", "w") as o1, \
                open('gene_predictions.gff3', 'w') as o2:
            for ev in evidence[0]:
                os.symlink(ev, os.path.basename(ev))
                o1.write('\t'.join(weights[os.path.basename(ev)]) + '\n')
                if os.path.basename(ev) \
                        not in ['nr_transcripts.gff3', 'prot.gff3']:
                    with open(ev, 'r') as i1:
                        o2.write(i1.read())
        # partition inputs to run in parallel
        _run(
            "partition_EVM_inputs.pl --genome " + self.args.genome +
            " --gene_predictions gene_predictions.gff3 " +
            "--transcript_alignments nr_transcripts.gff3 " +
            "--protein_alignments prot.gff3 %s " % evm_opts +
            "--segmentSize 100000 --overlapSize 10000 " +
            "--partition_listing partitions_list.out")
        # write a list of commands (to run in parallel on the splits)
        stops = '--stop_codons %s' % \
            ','.join(
                CodonTable.unambiguous_dna_by_id[int(self.args.code)].
                stop_codons)
        _run(
            "write_EVM_commands.pl --genome " + self.args.genome +
            " --weights %s " % os.path.join(os.getcwd(), 'weights.txt') +
            "--gene_predictions gene_predictions.gff3 " +
            "--protein_alignments prot.gff3 %s " % evm_opts +
            "--transcript_alignments nr_transcripts.gff3 " +
            "--output_file_name consensus_models.out %s " % stops +
            '--partitions partitions_list.out ', 'commands.list')
        cmds = [line.rstrip() for line in open('commands.list')]
        cmds_chunks = [
            cmds[x:x+int(self.numcpu)]
            for x in range(0, len(cmds), int(self.numcpu))]
        for chunk in cmds_chunks:
            procs = [
                subprocess.Popen('nice -19 %s' % cmd, shell=True)
                for cmd in chunk]
            for proc in procs:
                proc.wait()
        # now recombine all the partitions
        _run(
            "recombine_EVM_partial_outputs.pl " +
            "--partitions partitions_list.out " +
            "--output_file_name consensus_models.out")
        # make a file conversion
        _run(
            "convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out" +
            " --output consensus_models.out --genome " + self.args.genome)
        # gather all the consensi
        cons_files = [
            f for f in Path(os.getcwd()).rglob('consensus_models.out.gff3')]
        with open('consensus_models.gff3', 'w') as outfile:
            for f in cons_files:
                with open(f) as infile:
                    for line in infile:
                        outfile.write(line)
        return

    def add_UTRs_from_PASAdb(self, pasadb):
        """
        Takes consensus models and adds UTRs + models alternative splicing
        using the transcriptome database (created and loaded during assembly)
        """
        # create the required pasa files
        pasadb = str(Path('../%s' % pasadb).resolve())
        aligncfg = open('alignAssembly.config', 'w')
        aligncfg.write('DATABASE=%s\n' % pasadb)
        aligncfg.write(
            'validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=95\n')
        aligncfg.write('validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=95\n')
        aligncfg.write('subcluster_builder.dbi:-m=50\n')
        aligncfg.close()
        annotcfg = open('annotCompare.config', 'w')
        annotcfg.write('DATABASE=%s\n' % pasadb)
        annotcfg.close()
        _run(
            "Load_Current_Gene_Annotations.dbi -c alignAssembly.config " +
            "-g %s -P consensus_models.gff3" % self.args.genome)
        codes = {
            '6': '--GENETIC_CODE Tetrahymena',
            '10': '--GENETIC_CODE Euplotes',
            '12': '--GENETIC_CODE Candida',
            '11': ''}
        _run(
            "Launch_PASA_pipeline.pl -c annotCompare.config " +
            "-A -g %s --CPU %s " %
            (self.args.genome, self.numcpu) +
            " -t %s/%s %s" %
            (self.wd, self.args.transcriptsFasta, codes[self.args.code]))
        return

    def build_consensus_models(self, *evidence):
        """
        Routine to build consensus models from predictions using EVM
        """
        # run either evm or evigan
        self.run_evm(evidence)
        if self.args.utrs:
            # add UTRs if transcriptome assembly is available
            self.add_UTRs_from_PASAdb(self.args.utrs)
        pasa_out = glob.glob("*gene_structures_post_PASA_updates.*.gff3")
        cons = pasa_out[0] if pasa_out else 'consensus_models.gff3'
        # format the gff3 to properly represented nested features
        consdb = gffutils.create_db(
            cons, ':memory:', merge_strategy='create_unique')
        if os.path.exists('../orf_finder/transcript_orfs.gff3'):
            orf_transcripts = gffutils.create_db("../orf_finder/transcript_orfs.gff3", ":memory:")
            missing_models = gffintersecter.find_nonoverlapping_genes(orf_transcripts, consdb)
            consdb_patched = [f for f in consdb.all_features()]
            consdb_patched.extend(missing_models)
            consdb = gffutils.create_db(
                consdb_patched, ':memory:', merge_strategy='merge', from_string = True)
        consdb.dialect['order'].append('locus_tag')
        consdb = gffutils.create_db(gffparser.fix_CDS_phases(consdb), ':memory:', merge_strategy='merge')
        consdb = gffparser.prettify_gff3(consdb, self.shortname)
        with open('final.gff3', 'w') as outfile:
            for line in consdb:
                outfile.write('%s\n' % str(line))
        if os.path.exists('%s/final.gff3' % self.wd):
            os.remove('%s/final.gff3' % self.wd)
        shutil.copy('final.gff3', '%s/final.gff3' % self.wd)
        os.chdir(self.wd)
        sys.stdout.write(' done.\n')
        return


def _run(cmd, out_file=None):
    """"
    Run commands in the shell
    """
    with open(out_file if out_file else 'cmd.out', 'a') as out, \
        open('cmd.err', 'a') as err, \
            open('cmd.log', 'a') as log:
        res = subprocess.run(
            cmd.split(), universal_newlines=True, stdout=out, stderr=err)
        log.write("%s\ncmd: %s\n" % (datetime.datetime.now(), cmd))
        if res.stdout:
            if out_file:
                out.write(res.stdout)
            else:
                out.write(
                    "%s\ncmd: %s\n%s\n" %
                    (datetime.datetime.now(), cmd, res.stdout))
        if res.stderr:
            err.write(
                "%s\ncmd: %s\n%s\n" %
                (datetime.datetime.now(), cmd, res.stderr))
        if res.returncode > 0:
            raise Exception(
                'command:\n%s\nfailed to run. error:\n' % cmd)
    return


def _init(routine, output):
    if os.path.exists('%s/%s' % (routine, output)):
        try:
            gffutils.helpers.sanitize_gff_db(
                gffutils.create_db(
                    '%s/%s' % (routine, output), ':memory:', verbose=False))
            print('%s already done! Skipping.' % routine)
            return(True)
        except Exception as e:
            if re.search('No lines parsed', str(e)):
                print('ERROR: %s is empty or not in GFF3.' % output)
                print('The %s routine appears broken.' % routine)
                sys.exit(1)
    elif os.path.exists(routine):
        shutil.rmtree(routine)
        os.makedirs(routine)
        os.chdir(routine)
    else:
        os.makedirs(routine)
        os.chdir(routine)


def is_gff(gff):
    if os.path.exists(gff):
        try:
            fh = gffutils.create_db(gff, ':memory:')
            for f in fh.all_features():
                if len(f.attributes) == 0 or f.start is None or f.end is None:
                    print("File %s appears to be malformed gff\n" % gff)
                    sys.exit(1)
        except Exception as e:
            print(e)
            sys.exit(1)
    else:
        print('ERROR: %s does not exist.' % gff)
        sys.exit(1)


def is_fasta(fasta):
    if os.path.exists(fasta):
        try:
            with open(fasta, "r") as handle:
                fh = SeqIO.parse(handle, "fasta")
                if not any(fh):
                    print("Could not parse %s, file possibly not a proper fasta\n" % fasta)
        except Exception as e:
            print("Could not open %s\n" % fasta)
            print(e)
            sys.exit(1)
    else:
        print('ERROR: %s does not exist.' % fasta)
        sys.exit(1)


def _rand_string(length):
    l = string.ascii_uppercase
    result = ''.join(random.choice(l) for i in range(length))
    return(result)
