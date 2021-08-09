# Dockerfile for MutantHuntWGS pipeline
# Coded by Jen Walker and edited by Mitch Ellison

FROM ubuntu:latest

WORKDIR /Main



# remove interactive dialogue with debian front end during automated build

ENV DEBIAN_FRONTEND=noninteractive



# We need to get the basic bash commands we need

RUN apt-get -qq update && \
	apt-get -qq -y install apt-utils && \
	apt-get -qq -y install wget unzip bzip2 gcc git libncurses5-dev zlib1g-dev make python g++ libtbb-dev pkg-config default-jre



# Install MutantHuntWGS

RUN git clone https://github.com/mae92/MutantHuntWGS.git

RUN echo 'export PATH="/Main/MutantHuntWGS/Code:$PATH"' >> ~/.bashrc

RUN chmod 777 /Main/MutantHuntWGS/Code/MutantHuntWGS.sh

RUN chmod 777 /Main/MutantHuntWGS/Code/MutantHuntWGS_v1.1.sh



# Install bowtie2 (currently master, might want to grab 2.2.9)

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip/download && \
    unzip ./download && \
    cd ./bowtie2-2.2.9 && \
#    make && \
    cd ../

RUN echo 'export PATH="/Main/bowtie2-2.2.9:$PATH"' >> ~/.bashrc



# Install samtools 1.3.1

RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar -jxf ./samtools-1.3.1.tar.bz2 && \
    cd samtools-1.3.1 && \
    ./configure && \
    make install && \
    cd ../



# Install bcftools 1.3.1

RUN wget https://github.com/samtools/bcftools/releases/download/1.3/bcftools-1.3.tar.bz2 && \
    tar -jxf ./bcftools-1.3.tar.bz2 && \
    cd bcftools-1.3 && \
    make && \
    cd ../

RUN echo 'export PATH="/Main/bcftools-1.3:$PATH"' >> ~/.bashrc



# Install vcftools 0.1.14

RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz && \
    tar -zxf ./vcftools-0.1.14.tar.gz && \
    cd vcftools-0.1.14 && \
    ./configure && \
    make && \
    make install && \
    cd ../



# Install snpEff 4.3p

RUN wget http://sourceforge.net/projects/snpeff/files/snpEff_v4_3p_core.zip && \
    unzip ./snpEff_v4_3p_core.zip

RUN echo 'export PATH="/Main/snpEff/scripts:$PATH"' >> ~/.bashrc

RUN /Main/snpEff/scripts/snpEff ann -v Saccharomyces_cerevisiae



# Install SIFT4G annotator

RUN wget https://github.com/pauline-ng/SIFT4G_Annotator/raw/master/SIFT4G_Annotator.jar

RUN wget https://sift.bii.a-star.edu.sg/sift4g/public/Saccharomyces_cerevisiae/EF4.74.zip && \
	unzip ./EF4.74.zip 

RUN echo 'export PATH="/Main/SIFT4G_Annotator.jar:$PATH"' >> ~/.bashrc


ADD . .

RUN apt-get -qq -y install emacs

