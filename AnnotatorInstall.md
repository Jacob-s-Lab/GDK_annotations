# ANNOVAR
Use Ubuntu 20.4.
```bash
apt-get update --fix-missing

DEBIAN_FRONTEND=noninteractive apt-get --no-install-recommends install -y gcc g++ make wget zip git git-lfs perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev

apt-get clean
```

Download annovar.latest.tar.gz (Academic certification is required).
```bash
cd /opt
tar -zxvf annovar.latest.tar.gz
cd annovar

./annotate_variation.pl -buildver hg38 -downdb refGene humandb/
./annotate_variation.pl --buildver hg38 --downdb seq humandb/hg38_seq
./retrieve_seq_from_fasta.pl humandb/hg38_refGene.txt -seqdir humandb/hg38_seq -format refGene -outfile humandb/hg38_refGeneMrna.fa

```

---

# SnpEff
Use Ubuntu 20.4.
Requires Java 12 or higher.
```bash
apt-get update --fix-missing

DEBIAN_FRONTEND=noninteractive apt-get --no-install-recommends install -y \
openjdk-21-jdk gcc g++ make wget unzip git git-lfs

apt-get clean
```
Download SnpEff.
Version 5.2 (2023-09-29)
```bash
cd cd /opt
# wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
# unzip snpEff_latest_core.zip
wget https://github.com/pcingola/SnpEff/archive/refs/tags/v5.2.zip
```

```bash
unzip v5.2.zip
cd snpEff

java -jar snpEff.jar download -v GRCh38.mane.1.2.refseq
# (https://snpeff.blob.core.windows.net/databases/v5_2/snpEff_v5_2_GRCh38.mane.1.2.refseq.zip)

java -jar snpEff.jar download -v GRCh38.mane.1.2.ensembl
# (https://snpeff.blob.core.windows.net/databases/v5_2/snpEff_v5_2_GRCh38.mane.1.2.ensembl.zip)

```

---

# VEP
Use Ubuntu 20.4.
```bash
apt update --fix-missing
DEBIAN_FRONTEND=noninteractive apt-get --no-install-recommends install -y bzip2 parallel curl git make g++ wget

DEBIAN_FRONTEND=noninteractive apt-get --no-install-recommends install -y build-essential cpanminus libmysqlclient-dev zlib1g-dev libbz2-dev liblzma-doc liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev libpng-dev libexpat1-dev libxml-libxml-perl libbio-perl-perl

apt-get clean

cpanm --notest HTTP::Tiny Archive::Zip Archive::Extract DBD::mysql DBD::SQLite DBI LWP::Simple List::MoreUtils

SAMTOOLS_VERSION=1.18

cd /tmp \
  && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 -O samtools.tar.bz2 \
  && tar -xjvf samtools.tar.bz2 \
  && cd samtools-${SAMTOOLS_VERSION} \
  && make \
  && make install
cd ..
rm -rf /tmp/*
```

Download VEP version 111.0.
```bash
cd /opt
VEP_VERSION=111.0

wget https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/${VEP_VERSION}.tar.gz
tar -zxvf ${VEP_VERSION}.tar.gz
mv ensembl-vep-release-${VEP_VERSION} ensembl-vep

cd /opt/ensembl-vep

perl INSTALL.pl --AUTO a --NO_UPDATE --NO_TEST

mkdir /opt/VepCache

perl INSTALL.pl --AUTO c --NO_UPDATE --NO_TEST \
    --ASSEMBLY GRCh38 \
    --CACHE_VERSION 111 \
    --SPECIES homo_sapiens \
    --CACHEDIR /opt/VepCache/cache/

perl INSTALL.pl --AUTO c --NO_UPDATE --NO_TEST \
    --ASSEMBLY GRCh38 \
    --CACHE_VERSION 111 \
    --SPECIES homo_sapiens_refseq \
    --CACHEDIR /opt/VepCache/cache/

cd /opt/VepCache
curl -O https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
bgzip Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

