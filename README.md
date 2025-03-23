# Compare Three Annotators For Genetic Interpretation

---

# ClinVar Data Collection and Curation
Can Download our Normalized VCF and Decompress for following Annotation
[ClinVar_PLP](ClinVar_GRCh38_BLB_Tx1_20240107.vcf.gz)
[ClinVar_BLB](ClinVar_GRCh38_BLB_Tx1_20240107.vcf.gz)

OR See [ClinVar VCF Normalize](ClinVarVcfNorm.md) to Download Raw ClinVar VCF and Curation by yourself

Curated ClinVar Tsv File Download
[ClinVar_PLP](clinvar_GRCh38_PLP_Tx1.tsv.gz)
[ClinVar_BLB](clinvar_GRCh38_BLB_Tx1.tsv.gz)

OR See [ClinVar DataSet Curation](ClinVarSamplesetCuration.md)

---

# Annotators Installing
See [Annotator Installing](AnnotatorInstall.md)

---

# Variant Annotation

## ANNOVAR

Make ANNOVAR `avinput`
```bash
cd /opt/annovar

INPUT_DIR=/your/clinvar/vcf/path/ClinVar_VCF
OUTPUT_DIR=/your/output/volume/ClinVar20240107/refseq_annovar

perl convert2annovar.pl --includeinfo -format vcf4 $INPUT_DIR/ClinVar_GRCh38_PLP_Tx1_20240107.vcf > $OUTPUT_DIR/GRCh38_PLP/ClinVar_GRCh38_PLP.avinput

perl convert2annovar.pl --includeinfo -format vcf4 $INPUT_DIR/ClinVar_GRCh38_BLB_Tx1_20240107.vcf > $OUTPUT_DIR/GRCh38_BLB/ClinVar_GRCh38_BLB.avinput
```
Run ANNOVAR with RefSeq DB
thread number can change
```bash
IN_avinput=/your/output/volume/ClinVar20240107/refseq_annovar/GRCh38_PLP/ClinVar_GRCh38_PLP.avinput
OUT_DIR=/your/output/volume/ClinVar20240107/refseq_annovar/GRCh38_PLP

perl /volume/wgsa/annovar/table_annovar.pl $IN_avinput \
/volume/wgsa/annovar/humandb/ --buildver hg38 --protocol refGene --operation g -polish \
-out $OUT_DIR/ClinVar_GRCh38_PLP \
--thread 1 \
--argument '--neargene 5000 --hgvs --transcript_function --separate' \
--dot2underline --nastring .
```
```bash
IN_avinput=/your/output/volume/ClinVar20240107/refseq_annovar/GRCh38_BLB/ClinVar_GRCh38_BLB.avinput
OUT_DIR=/your/output/volume/ClinVar20240107/refseq_annovar/GRCh38_BLB

perl /volume/wgsa/annovar/table_annovar.pl $IN_avinput \
/volume/wgsa/annovar/humandb/ --buildver hg38 --protocol refGene --operation g -polish \
-out $OUT_DIR/ClinVar_GRCh38_BLB \
--thread 1 \
--argument '--neargene 5000 --hgvs --transcript_function --separate' \
--dot2underline --nastring .
```

## SnpEff
Run SnpEff with RefSeq DB
Xmx Ram can change
```bash
cd /opt/snpEff

INPUT_DIR=/your/clinvar/vcf/path/ClinVar_VCF
OUTPUT_DIR=/your/output/volume/ClinVar20240107/refseq_snpeff

SnpEff_INPUT_VCF=$INPUT_DIR/ClinVar_GRCh38_PLP_Tx1_20240107.vcf
SnpEff_OUTPUT_VCF=$OUTPUT_DIR/GRCh38_PLP/ClinVar_GRCh38_PLP_Tx1_20240107.snpeff.vcf

java -Xmx10g -jar snpEff.jar GRCh38.mane.1.2.refseq $SnpEff_INPUT_VCF > $SnpEff_OUTPUT_VCF
```
```bash
cd /opt/snpEff

INPUT_DIR=/your/clinvar/vcf/path/ClinVar_VCF
OUTPUT_DIR=/your/output/volume/ClinVar20240107/refseq_snpeff

SnpEff_INPUT_VCF=$INPUT_DIR/ClinVar_GRCh38_BLB_Tx1_20240107.vcf
SnpEff_OUTPUT_VCF=$OUTPUT_DIR/GRCh38_BLB/ClinVar_GRCh38_BLB_Tx1_20240107.snpeff.vcf

java -Xmx10g -jar snpEff.jar GRCh38.mane.1.2.refseq $SnpEff_INPUT_VCF > $SnpEff_OUTPUT_VCF
```

## VEP
Run VEP with RefSeq DB
BUFFER_SIZE, FORK can change
```bash
cd /opt/ensembl-vep

INPUT_DIR=/your/clinvar/vcf/path/ClinVar_VCF
OUTPUT_DIR=//your/output/volume/ClinVar20240107/refseq_vep

VEP_INPUT_VCF=$INPUT_DIR/ClinVar_GRCh38_PLP_Tx1_20240107.vcf
VEP_OUTPUT_VCF=$OUTPUT_DIR/GRCh38_PLP/ClinVar_GRCh38_PLP_Tx1_20240107.vep.vcf

ASSEMBLY=GRCh38
BUFFER_SIZE=500
CACHE_VERSION=111
FORK=2

./vep \
  --format vcf \
  -i $VEP_INPUT_VCF \
  --vcf \
  -o $VEP_OUTPUT_VCF \
  --assembly ${ASSEMBLY} \
  --refseq \
  --offline \
  --cache \
  --cache_version ${CACHE_VERSION} \
  --dir_cache /opt/VepCache \
  --fasta /opt/VepCache/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  --hgvsg \
  --hgvs \
  --total_length \
  --variant_class \
  --symbol \
  --protein \
  --canonical \
  --mane \
  --biotype \
  --no_stats \
  --no_escape \
  --buffer_size ${BUFFER_SIZE} \
  --fork ${FORK}
```
```bash
cd /opt/ensembl-vep

INPUT_DIR=/your/clinvar/vcf/path/ClinVar_VCF
OUTPUT_DIR=//your/output/volume/ClinVar20240107/refseq_vep

VEP_INPUT_VCF=$INPUT_DIR/ClinVar_GRCh38_BLB_Tx1_20240107.vcf
VEP_OUTPUT_VCF=$OUTPUT_DIR/GRCh38_PLP/ClinVar_GRCh38_BLB_Tx1_20240107.vep.vcf

ASSEMBLY=GRCh38
BUFFER_SIZE=500
CACHE_VERSION=111
FORK=2

./vep \
  --format vcf \
  -i $VEP_INPUT_VCF \
  --vcf \
  -o $VEP_OUTPUT_VCF \
  --assembly ${ASSEMBLY} \
  --refseq \
  --offline \
  --cache \
  --cache_version ${CACHE_VERSION} \
  --dir_cache /opt/VepCache \
  --fasta /opt/VepCache/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  --hgvsg \
  --hgvs \
  --total_length \
  --variant_class \
  --symbol \
  --protein \
  --canonical \
  --mane \
  --biotype \
  --no_stats \
  --no_escape \
  --buffer_size ${BUFFER_SIZE} \
  --fork ${FORK}
```

---

