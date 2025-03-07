#!/bin/bash
##############################################################
##############################################################
# PREAMBLE
##############################################################
##############################################################
# stop function
stop () {
	echo -n "Pipeline paused, continue with <return>"
	read REPLY
}

usage () { 
	echo
	printf "Usage: %s [OPTIONS] <ARGS>\n
	[OPTIONS] and corresponding <ARGS> are:\n
	Either paired-end:
		[-l] <left reads>
		[-r] <right reads>
	or single-end:
		[-s] <single-end reads>
	[-m] <min intron length> # default 20
	[-M] <max intron length> # default 5000
	[-g] <genome fasta>
	[-p] <phred quality score (33 for MISEQ, 64 for HISEQ)> # default 33
	[-n] <number of CPUs> # default is MAX
	[-S] <specificy strand-specific assembly, either RF or FR> # default off, i.e. unstranded
	[-A] <switch on read mapping>
	[-E] <switch to extract reads>
	[-T] <switch on Trinity assembly>
	[-e] <switch on StringTie assembly>
	[-P] <switch on PASA alignment>
	[-c] <genetic code according to ncbi table>
	[-h] Display this help message
	[-j] switch on jaccard clipping (for gene-dense orgnanisms and high coverage data)
	[-t] <EndToEnd/Local> # default Local
	\n", "$0" 1>&2; exit 1; 
}
# Make sure there are arguments supplied to the script
if [ $# -lt 1 ] ; then
	echo No options were supplied to the script. See usage.
		usage
fi
map_reads () { # STAR read mapping
	echo Mapping reads to the genome
	echo 1: building genome index
	if [ ! -d build-index ]; then
		mkdir build-index
		nice -19 STAR --genomeSAindexNbases 3 --limitGenomeGenerateRAM 40317074816 --runThreadN ${numthreads} --runMode genomeGenerate --genomeDir build-index --genomeFastaFiles ${genome}
		# eval nice -19 STAR --limitGenomeGenerateRAM 40317074816 --runThreadN "$numthreads" --runMode genomeGenerate --genomeDir build-index --genomeFastaFiles "$genome"
	fi
	echo 2: mapping reads reads
	[ "$maxIntronLen" != "" ] && STARmaxIntronLen="--alignIntronMax $maxIntronLen" || STARmaxIntronLen=
	filetype=$(echo $reads | sed 's/[, ]/\n/g' | xargs file -L | awk '{print $2}' | sort | uniq)
	[ "$filetype" = "gzip" ] && zcat="--readFilesCommand zcat"
	nice -19 STAR --runThreadN ${numthreads} --genomeDir build-index ${alignEndsType} --readFilesIn ${reads} --outSAMtype BAM SortedByCoordinate --outSJfilterIntronMaxVsReadN 100 300 500 --alignIntronMin ${minIntronLen} ${STARmaxIntronLen} --outFileNamePrefix STAR_ --outSAMattributes All --outSAMattrIHstart 0 --outSAMstrandField intronMotif --limitBAMsortRAM 27643756136 ${zcat} ${quality}
	# eval nice -19 STAR --runThreadN "$numthreads" --genomeDir build-index "$alignEndsType" --readFilesIn "$reads" --outSAMtype BAM SortedByCoordinate --alignIntronMin "$minIntronLen" "$STARmaxIntronLen" --outFileNamePrefix STAR_ --outSAMattributes Standard --limitBAMsortRAM 27643756136 "$zcat" "$quality"
	starout="$?"
	if [ "$starout" -ne 0 ]; then
		echo "STAR command failed, maybe read lengths are not compatible with that version? Trying STARlong"
		nice -19 STARlong --runThreadN ${numthreads} --genomeDir build-index ${alignEndsType} --readFilesIn ${reads} --outSAMtype BAM SortedByCoordinate --outSJfilterIntronMaxVsReadN 100 300 500 --alignIntronMin ${minIntronLen} ${STARmaxIntronLen} --outFileNamePrefix STAR_ --outSAMattributes All --outSAMattrIHstart 0 --outSAMstrandField intronMotif --limitBAMsortRAM 11391164255 ${zcat} ${quality}
	fi
	awk 'BEGIN{FS=OFS="\t"} {if ($4==2) {$4="-"} else if ($4==1) {$4="+"} else {$4="."}; if ($5==0) {$5="non-canonical"} else if ($5==1) {$5="GT_AG"} else if ($5==2) {$5="CT_AC"} else if ($5==3) {$5="GC_AG"} else if ($5==4) {$5="CT_GC"} else if ($5==5) {$5="AT_AC"} else {$5="GT_AT"}; print $1,"STAR","intron",$2,$3,$7+$8,$4,".","mult="$7+$8";pri=4;src=E"}' STAR_SJ.out.tab > hints_introns.gff;
	samtools view -b -f 0x10 STAR_Aligned.sortedByCoord.out.bam > STAR_reverse.bam;
	samtools view -b -F 0x10 STAR_Aligned.sortedByCoord.out.bam > STAR_forward.bam;
	bam2wig STAR_reverse.bam > minus.wig;
	bam2wig STAR_forward.bam > plus.wig;
	wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=exonpart --radius=4.5 --pri=4 --strand="-" < minus.wig > hints.ep.minus.gff;
	wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=exonpart --radius=4.5 --pri=4 --strand="+" < plus.wig > hints.ep.plus.gff;
	cat hints.ep.minus.gff hints.ep.plus.gff > hints_coverage.gff;
	rm -rf build-index STAR_reverse.bam STAR_forward.bam minus.wig plus.wig hints.ep.minus.gff hints.ep.plus.gff
	echo done.
}
trinity_assembly () { # run the genome-guided trinity assembly
	echo Running assembly machines...
	echo 1: running genome-guided Trinity
	if [ ! -f trinity-gg.fasta ]; then
		nice -19 Trinity --genome_guided_bam STAR_Aligned.sortedByCoord.out.bam --genome_guided_max_intron ${maxIntronLen} --max_memory ${mem} --CPU ${numthreads} --full_cleanup --output trinity-gg --SS_lib_type ${libtypeTrinity} ${jaccard}
		mv trinity-gg/Trinity-GG.fasta trinity-gg.fasta && rm -rf trinity-gg
		rm -rf STAR_Aligned.sortedByCoord.out.bam.*
	else
		echo Genome-guided assembly already done, skipping
	fi
	echo 2: running denovo Trinity
	if [ ! -f trinity-denovo.fasta ]; then
		nice -19 Trinity --seqType fq --max_memory ${mem} ${reads} --CPU ${numthreads} --output trinity-denovo --full_cleanup --SS_lib_type ${libtypeTrinity} ${jaccard}
		mv trinity-denovo-mapped-reads.Trinity.fasta trinity-denovo-mapped-reads.fasta && rm -rf trinity-denovo-mapped-reads
		mv trinity-denovo.Trinity.fasta trinity-denovo.fasta && rm -rf trinity-denovo
	else
		echo De novo assembly already done, skipping
	fi
	echo done.
}
run_PASA () { # run PASA on the assemblies
	echo Assembling spliced alignments from the transcriptome assemblies
	rm -rf alignment.validations.output ./*cidx blat* gmap* ./*gmap 11.ooc pasa_run.* cleaning_* trinity-comprehensive* outparts_* seqcl* err_seqcl* tdn.accs cleaning_* 
	# create the pasa configs (to run PASA)
	db=$(pwd)/$name.sqlite
	awk -v db=${db} '$1=="DATABASE=" {$1="DATABASE="db}1' ${PASAHOME}/annotCompare.config > annotCompare.config
	awk -v db=${db} '$1=="DATABASE=" {$1="DATABASE="db}1' ${PASAHOME}/alignAssembly.config > alignAssembly.config
	# pull out accessions
	awk '/^>/ {gsub(/>/, "", $1); print $1}' trinity-denovo.fasta > tdn.accs
	# cat the output from the GG/denovo assemblies
	cat trinity-denovo.fasta trinity-gg.fasta > trinity-comprehensive.fasta
	# clean the DB
	echo 1: trim poly-A tails
	nice -19 seqclean trinity-comprehensive.fasta -l 90 -c 6
	echo 2: running program to assemble spliced alignments on the comprehensive transcriptome assembly
	Launch_PASA_pipeline.pl -c alignAssembly.config -C -r -R -g ${genome} -t trinity-comprehensive.fasta.clean -T -u trinity-comprehensive.fasta --ALIGNERS gmap,blat --CPU ${numthreads} --TDN tdn.accs -I ${maxIntronLen} --stringent_alignment_overlap 30.0 ${strand} ${gencode} ${pasastringtie}
	echo 3: build a non-redundant comprehensive transcriptome assembly
	build_comprehensive_transcriptome.dbi -c alignAssembly.config -t ${name}.sqlite.assemblies.fasta --min_per_ID 95 --min_per_aligned 95 
	# remove redundancy in PASA fasta file
	cp compreh_init_build/compreh_init_build.fasta compreh_init_build/compreh_init_build.gff3 .
	awk 'BEGIN{RS=">"; ORS=">"} !x[$0]++' compreh_init_build.fasta | sed '$d' > nr_transcripts.fasta
	awk 'BEGIN{FS=OFS="\t"}{split($9,a,"[ ;=]"); ++count;$2="PASA-assembly";$3="exon";$9="ID="a[4]":exon:"count";Parent="a[4]";";print}' compreh_init_build.gff3 | tee nr_transcripts.gff3 | awk 'BEGIN{FS=OFS="\t"}{split($9,a,"[=;]"); $9="pri=3;src=E;group="a[4]}1' > hints_transcripts.gff
	cat hints_transcripts.gff hints_introns.gff hints_coverage.gff > hints_rnaseq.gff
	#rm -rf alignment.validations.output ./*cidx ./*blat* ./*gmap* 11.ooc pasa_run.* cleaning_* trinity-comprehensive* outparts_* seqcl* err_seqcl* tdn.accs ./*.err Log.out
	echo done.
}

##############################################################
##############################################################
# MAIN
##############################################################
##############################################################
# set some reasonable defaults
minIntronLen=20
jaccard=
gencode="--GENETIC_CODE universal"
mem=$(free -g | grep '^Mem:' | awk '{printf "%0.fG",$2/2}')
numthreads=$(lscpu | awk '/^CPU\(s\):/ {print $2/2}') # default number of CPUs
zcat=
alignEndsType="--alignEndsType Local"

# parse the options from the command line
while getopts 'n:g:l:r:es:m:M:p:S:hAEPt:Tj' OPTION
do
	case $OPTION in
		g)	fc "$OPTARG"
			genome="$OPTARG"
			name=$(echo "$genome" | cut -f1 -d'.');;
		l)	fc "$OPTARG"
			left="$OPTARG";;
		r)	fc "$OPTARG"
			right="$OPTARG";;
		s)	fc "$OPTARG"
			singles="$OPTARG";;
		m)	minIntronLen="$OPTARG";;
		M)	maxIntronLen="$OPTARG";;
		p)	[ "$OPTARG" -eq 64 ] && quality="--outQSconversionAdd -31";;
		S)	libtypeTrinity="$OPTARG"; libtypeStringTie=$(echo "$OPTARG" | awk '{print "--"tolower($1)}');;
		t)	alignEndsType="--alignEndsType $OPTARG";;
		A)	readmapping="true";;
		E)	extractreads="true";;
		T)	runtrinity="true";;
		P)	runpasa="true";;
		n)	numthreads="$OPTARG";;
		j)	jaccard=" --jaccard_clip";;
		j)	[ "$OPTARG" -eq 6 ] && gencode=" --GENETIC_CODE Tetrahymena";;
		h)	usage;;
		?)	usage;;
	esac
done
shift $((OPTIND - 1))

if [ -n "$readmapping" ]; then
	[ "$left" != "" -a "$right" != "" ] && reads="$left $right" || reads="$singles"
	map_reads
fi
#maxIntronLen=$(awk '{print $5-$4}' hints_introns.gff | Rscript -e 'tmp <- as.numeric (readLines ("stdin")); median(tmp)+round(3*sd(tmp),0)' | awk '{print $2}')
if [ -n "$extractreads" ]; then
	mapped_read_extraction
fi
if [ -n "$runtrinity" ]; then
	[ "$left" != "" -a "$right" != "" ] && reads="--left $left --right $right" || reads="--single $singles"
	trinity_assembly
fi
if [ -n "$runpasa" ]; then
	[ -n "$libtypeTrinity" ] && strand="--transcribed_is_aligned_orient" || strand=""
	run_PASA
fi
if [ -z "$readmapping" -a -z "$extractreads" -a -z "$runtrinity" -a -z "$runpasa" ]; then
	map_reads
	mapped_read_extraction
	trinity_assembly
	run_PASA
fi

exit 0
