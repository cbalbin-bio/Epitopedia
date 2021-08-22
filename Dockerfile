FROM ubuntu:20.04
ENV LISTEN_PORT=5000
EXPOSE 5000
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && apt-get install -y \
    build-essential \
    autotools-dev \
    automake \
    autoconf \
    libbz2-dev \
    libboost-all-dev \
    python3.9 \
    python3-pip \
    sqlite3 \
    wget \
    curl \
    gzip \
    libssl-dev \
    libio-socket-ssl-perl \
    cpanminus \
    libxml-simple-perl \ 
    libwww-perl \
    libnet-perl \
    golang


RUN python3.9 -m pip install flask gemmi biopython rich dataclasses-json


WORKDIR /app

# Get blast
RUN wget "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz" ;\
    tar -zxvf ncbi-blast-2.12.0+-x64-linux.tar.gz



# build DSSP
RUN wget "https://github.com/cmbi/dssp/archive/refs/tags/2.3.0.tar.gz" ;\
    tar -zxvf 2.3.0.tar.gz
WORKDIR /app/dssp-2.3.0
RUN ./autogen.sh ;\
    ./configure ;\
    make ;\
    make install
WORKDIR /app
RUN rm -r dssp-2.3.0 2.3.0.tar.gz

# build TMalign
RUN wget "https://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.f" ;\
    gfortran -static -O3 -ffast-math -lm -o TMalign TMalign.f ;\
    mv TMalign /usr/local/bin ;\
    rm TMalign.f

# get mysql2sqlite converter to convert IEDB's mySQL DB
RUN wget "https://raw.githubusercontent.com/dumblob/mysql2sqlite/d14d22ad7029cdf4d11825ee3c96922e8fbb0122/mysql2sqlite" ;\
    chmod +x mysql2sqlite ;\
    mv mysql2sqlite /usr/local/bin

# get mmseqs2
RUN wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz ;\
    tar xvfz mmseqs-linux-sse41.tar.gz ;\
    rm mmseqs-linux-sse41.tar.gz






RUN curl -s "ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz" -o "edirect.tar.gz" ;\
    tar -xf edirect.tar.gz ;\
    rm edirect.tar.gz

WORKDIR /app/edirect

RUN wget ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/xtract.Linux.gz ;\
    gunzip -f xtract.Linux.gz ;\
    chmod +x xtract.Linux;\
    wget ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/transmute.Linux.gz ;\
    gunzip -f transmute.Linux.gz ;\
    chmod +x transmute.Linux ;\
    wget ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/rchive.Linux.gz ;\
    gunzip -f rchive.Linux.gz ;\
    chmod +x rchive.Linux 
# RUN useradd -rm -d /home/ubuntu -s /bin/bash -G sudo -u 1001 ubuntu
# USER ubuntu
# add mmseqs2 to path

WORKDIR /app
COPY . /app

# ENV PATH=/usr/local/ncbi/edirect:${PATH}

RUN chmod +x run_epitopedia.py ; chmod +x generate_database.py

ENV PATH=/app/edirect:/app/mmseqs/bin/:/app/ncbi-blast-2.12.0+/bin:/app/:$PATH





