#!/usr/bin/env bash

set -exvuo pipefail

SENTIEON_VERSION=202308.01
SENTIEON_URL="https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-"
SENTIEON_ARM_URL="https://s3.amazonaws.com/sentieon-release/software/arm-sentieon-genomics-"

INSTALL_DIR="/home/azureuser/programs"
WORK_DIR="/home/azureuser/work"
mkdir -p "$INSTALL_DIR" "$WORK_DIR"


# Update the machine and install some packages
sudo apt update -y
sudo apt install -y
sudo apt install -y autoconf automake make gcc libdata-dump-perl zlib1g zlib1g-dev  \
    bzip2 libbz2-dev xz-utils curl openssl libncurses5-dev libncursesw5-dev libtool \
    nasm yasm gzip git unzip hostname cmake make gcc g++ autoconf liblzma-dev

# Get the software packge
cd "$INSTALL_DIR"
mkdir release
cd release
architechure=$(uname -m)
if [[ "$architechure" == "x86_64" ]]; then
    curl -L "${SENTIEON_URL}${SENTIEON_VERSION}.tar.gz" | tar -zxf -
    sentieon-genomics-"$SENTIEON_VERSION"/bin/sentieon driver -h
elif [[ "$architechure" == "arm" || "$architechure" == "aarch64" ]]; then
    curl -L "${SENTIEON_ARM_URL}${SENTIEON_VERSION}.tar.gz" | tar -zxf -
    arm-sentieon-genomics-"$SENTIEON_VERSION"/bin/sentieon driver -h
else
    echo "UNKNOWN ARCHITECHURE"
    exit 1
fi

# install jemalloc
cd "$INSTALL_DIR"
mkdir -p jemalloc
cd jemalloc
curl -L https://github.com/jemalloc/jemalloc/releases/download/5.2.1/jemalloc-5.2.1.tar.bz2 | \
    tar -jxf -
cd jemalloc-5.2.1/
./configure
make

# install igzip
cd "$INSTALL_DIR"
mkdir -p isa-l
cd isa-l
curl -L https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz | \
    tar -zxf -
cd isa-l-2.30.0/
./autogen.sh
./configure --prefix=/usr --libdir=/usr/lib
make && sudo make install

# install seqtk
cd "$INSTALL_DIR"
mkdir -p seqtk
cd seqtk
git clone https://github.com/lh3/seqtk.git
cd seqtk
make

# install pigz
cd "$INSTALL_DIR"
mkdir -p pigz
cd pigz
curl -L https://zlib.net/pigz/pigz-2.8.tar.gz | tar -zxf -
cd pigz-2.8
make

# install bedtools
cd "$INSTALL_DIR"
mkdir -p bedtools/2.30.0
cd bedtools/2.30.0
if [[ "$architechure" == "x86_64" ]]; then
    curl -L -o "bedtools" "https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary"
    chmod ugo+x ./bedtools
elif [[ "$architechure" == "arm" || "$architechure" == "aarch64" ]]; then
     wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
     tar -zxvf bedtools-2.30.0.tar.gz
     cd bedtools2
     sed -i '/^#include <cstdlib>/s/^/#include <cstdint>\n/' src/utils/general/ParseTools.h
     make
     mv bin/bedtools ../2.30.0/
     
else
    echo "UNKNOWN ARCHITECHURE"
    exit 1
fi

# install samtools
cd "$INSTALL_DIR"
mkdir -p samtools
cd samtools
curl -L "https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2" | \
    tar -jxf - 
cd samtools-1.15.1
./configure && make

# install bcftools
cd "$INSTALL_DIR"
mkdir -p bcftools
cd bcftools
curl -L "https://github.com/samtools/bcftools/releases/download/1.15.1/bcftools-1.15.1.tar.bz2" | \
    tar -jxf -
cd bcftools-1.15.1
./configure && make

# install rtgtools
cd "$INSTALL_DIR"
mkdir -p rtgtools
cd rtgtools
curl -L -O "https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.12.1/rtg-tools-3.12.1-nojre.zip"
unzip rtg-tools-3.12.1-linux-x64.zip
## Configure rtgtools
rtg_config="rtg-tools-3.12.1/rtg.cfg"
echo "RTG_TALKBACK=" > "$rtg_config"
echo "RTG_USAGE=" >> "$rtg_config"
echo "RTG_USAGE_OPTIONAL=username,hostname" >> "$rtg_config"

# Install snakemake with pip
cd "$INSTALL_DIR"
if [[ "$architechure" == "x86_64" ]]; then
    curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b
elif [[ "$architechure" == "arm" || "$architechure" == "aarch64" ]]; then
    curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh
    bash Miniconda3-latest-Linux-aarch64.sh -b
else
    echo "UNKNOWN ARCHITECHURE"
    exit 1
fi

python3=/home/azureuser/miniconda3/bin/python
pip3=/home/azureuser/miniconda3/bin/pip
$pip3 install snakemake pandas
snakemake=/home/azureuser/miniconda3/bin/snakemake

# Install the sentieon-cli
$pip3 install poetry
poetry=/home/azureuser/miniconda3/bin/poetry
cd "$INSTALL_DIR"
git clone https://github.com/Sentieon/sentieon-cli.git
cd sentieon-cli
$poetry config virtualenvs.in-project true
$poetry install

# install minimap2 for arm
curl -L "https://github.com/lh3/minimap2/archive/refs/tags/v2.26.tar.gz" | \
    tar -zxf -
cd minimap2-2.26
make arm_neon=1 aarch64=1

# install miniconda2 for hap.py
cd "$INSTALL_DIR"
if [[ "$architechure" == "x86_64" ]]; then
    curl -LO https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
    bash Miniconda2-latest-Linux-x86_64.sh -b
    python2=/home/azureuser/miniconda2/bin/python
    pip2=/home/azureuser/miniconda2/bin/pip
    pip2 install pysam pandas scipy
    # install happy
    cd "$INSTALL_DIR"
    mkdir -p happy
    cd happy
    curl -L "https://github.com/Illumina/hap.py/archive/refs/tags/v0.3.15.tar.gz" | \
        tar -zxf -
    cd hap.py-0.3.15
    python2 ./install.py "$INSTALL_DIR"/happy/hap.py-install --no-tests
fi