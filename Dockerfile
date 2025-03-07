# adapated parts from https://github.com/chrishah/maker-docker and https://github.com/Gaius-Augustus/Augustus/blob/master/Dockerfile

FROM ubuntu:18.04 AS annot-base

# install required packages
RUN apt-get --fix-missing update
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
RUN apt-get install -y build-essential wget autoconf unzip language-pack-en git default-jre \
  cmake \
	libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev \
	libsqlite3-dev libmysql++-dev \
	libboost-iostreams-dev zlib1g-dev \
	libbamtools-dev \
	libbz2-dev liblzma-dev \
	libncurses5-dev apt-utils \
	libssl-dev libcurl3-dev \
	python3-biopython \
	perl bioperl python3 python3-pip \
	exonerate ncbi-blast+ \
	ncbi-blast+-legacy cdbfasta \
	python-biopython python-pip python-gtk2 liblbfgs-dev libgsl-dev \
	parallel libopenmpi-dev openmpi-bin openmpi-common \
	&& ln -s /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgsl.so.0

WORKDIR /root

# build htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && tar jxvf htslib-1.11.tar.bz2 -C /root/ && mv /root/htslib-1.11 /root/htslib
WORKDIR "/root/htslib"
RUN autoheader && autoconf && ./configure && make -j2 && make install
# build bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 && tar jxvf bcftools-1.11.tar.bz2 -C /root/ && mv /root/bcftools-1.11 /root/bcftools
WORKDIR "/root/bcftools"
RUN autoheader && autoconf && ./configure && make -j2 && make install
# build samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && tar jxvf samtools-1.11.tar.bz2 -C /root/ && mv /root/samtools-1.11 /root/samtools
WORKDIR "/root/samtools"
RUN autoheader && autoconf -Wno-syntax && ./configure && make -j2 && make install
ENV TOOLDIR="/root"

# build augustus
WORKDIR "/root"
RUN git clone https://github.com/Gaius-Augustus/Augustus
RUN cd Augustus && git checkout 0e2e3114b0cade36e9b68398f1cdcb6bf5bdabe1 && sed -i '/TOOLDIR=/ s/\$(HOME)/root/' /root/Augustus/auxprogs/bam2wig/Makefile
COPY eukka.MPEW.RM.cfg /root/Augustus/config/extrinsic/
WORKDIR "/root/Augustus/"
RUN make -j2 && make install
ENV PATH="/root/Augustus/bin:${PATH}"
ENV PATH="/root/Augustus/scripts:${PATH}"
ENV AUGUSTUS_CONFIG_PATH /root/Augustus/config

# build snap
WORKDIR "/root"
RUN git clone https://github.com/KorfLab/SNAP
WORKDIR "/root/SNAP"
RUN make -j2
ENV PATH="/root/SNAP:${PATH}"
ENV ZOE="/root/SNAP/Zoe"

#install perl modules required for Genemark
RUN cpan YAML && \
        cpan Hash::Merge && \
        cpan Logger::Simple && \
        cpan Parallel::ForkManager && mkdir /root/Genemark
#add the (at this point empty) Genemark directory to the path in anticipation
ENV PATH="/root/Genemark:${PATH}"

# build codingquarry
RUN wget https://master.dl.sourceforge.net/project/codingquarry/CodingQuarry_v2.0.tar.gz && tar zxvf CodingQuarry_v2.0.tar.gz -C /root/
WORKDIR "/root/CodingQuarry_v2.0"
RUN make -j2
ENV PATH="/root/CodingQuarry_v2.0:${PATH}"
ENV QUARRY_PATH="/root/CodingQuarry_v2.0/QuarryFiles"

# add blast executables to location expected by repeatmodeler
RUN for f in $(find /usr/bin/ -name '*blast*'); do ln -s $f /usr/local/bin/; done
RUN mkdir /root/trf && cd /root/trf
RUN wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 && mv trf*.linux64 trf && chmod +x trf
ENV PATH="/root/trf:${PATH}"

# install genemark
COPY gmes_linux_64.tar.gz /root
COPY gm_key_64.gz /root
RUN tar zxvf /root/gmes_linux_64.tar.gz -C /root
RUN gunzip /root/gm_key_64.gz -c > ~/.gm_key
ENV PATH="/root/gmes_linux_64:${PATH}"

# install hmmer
RUN wget http://eddylab.org/software/hmmer/hmmer-3.3.1.tar.gz && tar zxvf hmmer-3.3.1.tar.gz -C /opt
WORKDIR "/opt/hmmer-3.3.1"
RUN ./configure && make -j2 && make install

# trinity (installed before gmap given conflict in except.h)
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.11.0/trinityrnaseq-v2.11.0.FULL.tar.gz && tar zxvf trinityrnaseq-v2.11.0.FULL.tar.gz -C /opt
WORKDIR "/opt/trinityrnaseq-v2.11.0"
RUN make -j2
RUN make plugins
ENV TRINITY_HOME="/opt/trinityrnaseq-v2.11.0"
ENV PATH="/opt/trinityrnaseq-v2.11.0:${PATH}"

# install blat for pasa
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat -O /usr/local/bin/blat && chmod +x /usr/local/bin/blat
# install evidencemodeler
RUN wget https://github.com/EVidenceModeler/EVidenceModeler/archive/v1.1.1.tar.gz && tar zxvf v1.1.1.tar.gz -C /root
ENV PERL5LIB="/root/EVidenceModeler-1.1.1/PerlLib:${PERL5LIB}"
ENV PATH="/root/EVidenceModeler-1.1.1/EvmUtils:${PATH}"
ENV PATH="/root/EVidenceModeler-1.1.1/EvmUtils/misc:${PATH}"
ENV EVM_HOME="/root/EVidenceModeler-1.1.1"

# jellyfish
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz && tar zxvf jellyfish-2.3.0.tar.gz -C /root
WORKDIR "/root/jellyfish-2.3.0/"
RUN ./configure && make && make install

# install gsnap/gmap for pasa
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2017-11-15.tar.gz && tar xzvf gmap-gsnap-2017-11-15.tar.gz -C /root
WORKDIR "/root/gmap-2017-11-15"
RUN ./configure && make -j2 && make install
# install fasta3 for pasa
RUN wget https://github.com/wrpearson/fasta36/releases/download/fasta-v36.3.8g/fasta-36.3.8g-linux64.tar.gz && tar zxvf fasta-36.3.8g-linux64.tar.gz -C /root
WORKDIR "/root/fasta-36.3.8g/src"
RUN make -f ../make/Makefile.linux_sse2 all && cp ../bin/fasta36 ../bin/fasta
ENV PATH="/root/fasta-36.3.8g/bin:${PATH}"

# install pasa
RUN apt-get install -y sqlite libdbd-sqlite3 libdbi-perl libgd-perl
#RUN cpan GD DBI
# pasapipeline v2.4.1
RUN wget https://github.com/PASApipeline/PASApipeline/releases/download/pasa-v2.4.1/PASApipeline.v2.4.1.FULL.tar.gz && tar zxvf PASApipeline.v2.4.1.FULL.tar.gz -C /root
WORKDIR "/root/PASApipeline.v2.4.1/"
RUN make -j2
ENV PASAHOME /root/PASApipeline.v2.4.1
ENV PATH="/root/PASApipeline.v2.4.1:${PATH}"
ENV PATH="/root/PASApipeline.v2.4.1/bin:${PATH}"
ENV PATH="/root/PASApipeline.v2.4.1/misc_utilities:${PATH}"
ENV PATH="/root/PASApipeline.v2.4.1/scripts:${PATH}"
COPY annotCompare.config /root/PASApipeline.v2.4.1/
COPY alignAssembly.config /root/PASApipeline.v2.4.1/

# STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz && tar zxvf 2.7.6a.tar.gz -C /opt
ENV PATH="/opt/STAR-2.7.6a/bin/Linux_x86_64_static:${PATH}"

# spaln
RUN wget https://github.com/ogotoh/spaln/archive/Ver.2.4.2.tar.gz && tar zxvf Ver.2.4.2.tar.gz -C /root && cd /root/spaln-Ver.2.4.2/src && ./configure && make && make install
ENV PATH="/root/spaln-Ver.2.4.2/bin:${PATH}"
ENV PATH="/root/spaln-Ver.2.4.2/seqdb:${PATH}"
ENV ALN_TAB="/root/spaln-Ver.2.4.2/table"
ENV ALN_DBS="/root/spaln-Ver.2.4.2/seqdb"

# fitild (compiled with g++-4.7 since code is old)
COPY fitild.tar.gz /root
RUN tar zxvf /root/fitild.tar.gz -C /root/ && cat /root/fitild/table/IldModel*.txt > /root/spaln-Ver.2.4.2/table/IldModel.txt
ENV PATH="/root/fitild/bin:${PATH}"
RUN chmod +x /root/PASApipeline.v2.4.1/bin/*

FROM annot-base

ARG USER_ID
ARG GROUP_ID

# fitild (compiled with g++-4.7 since code is old)
COPY fitild.tar.gz /opt
RUN tar zxvf /opt/fitild.tar.gz -C /opt/ && cat /opt/fitild/table/IldModel*.txt > /opt/spaln-Ver.2.4.2/table/IldModel.txt
ENV PATH="/opt/fitild/bin:${PATH}"

COPY . /opt/eukan
ENV PATH="/opt/eukan:${PATH}"
RUN python3 -m pip install -r /opt/eukan/requirements.txt

RUN chown -R $USER_ID:$GROUP_ID /opt/spaln-Ver.2.4.2 /opt/fitild

SHELL ["/bin/bash", "-c"]
