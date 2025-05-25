# Compare Three Annotators For Genetic Interpretation

---

# ClinVar Data Collection and Curation
You can download our normalized VCF file and decompress it for annotation:
- [ClinVar_PLP](ClinVar_GRCh38_BLB_Tx1_20240107.vcf.gz)
- [ClinVar_BLB](ClinVar_GRCh38_BLB_Tx1_20240107.vcf.gz)

Or, See [ClinVar VCF Normalize](ClinVarVcfNorm.md) for instructions on how to download the raw ClinVar VCF and perform curation yourself.

Curated ClinVar TSV files are also available for download:
- [ClinVar_PLP](clinvar_GRCh38_PLP_Tx1.tsv.gz)
- [ClinVar_BLB](clinvar_GRCh38_BLB_Tx1.tsv.gz)

Alternatively, see [ClinVar DataSet Curation](ClinVarSamplesetCuration.md) for more details.

---

# Annotators Installing
See [Annotator Installing](AnnotatorInstall.md) for installation instructions.

---

# Variant Annotation

## ANNOVAR

Prepare ANNOVAR `avinput` files:
```bash
cd /opt/annovar

INPUT_DIR=/your/clinvar/vcf/path/ClinVar_VCF
OUTPUT_DIR=/your/output/volume/ClinVar20240107/refseq_annovar

perl convert2annovar.pl --includeinfo -format vcf4 $INPUT_DIR/ClinVar_GRCh38_PLP_Tx1_20240107.vcf > $OUTPUT_DIR/GRCh38_PLP/ClinVar_GRCh38_PLP.avinput

perl convert2annovar.pl --includeinfo -format vcf4 $INPUT_DIR/ClinVar_GRCh38_BLB_Tx1_20240107.vcf > $OUTPUT_DIR/GRCh38_BLB/ClinVar_GRCh38_BLB.avinput
```
Run ANNOVAR with the RefSeq database (you can adjust the number of threads as needed):
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
Run SnpEff with the RefSeq database (you can adjust the Xmx RAM value as needed):
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
Run VEP with the RefSeq database (BUFFER_SIZE and FORK can be adjusted as needed):
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

# Trrascript Curation After Annotation
We provide an R script to process and split transcripts following annotation by annotators. You can then use the supplied `key` together with the processed `RefSeq transcript (NM)` identifiers as join criteria for comparative analyses in the programming language of your choice.
- [ANNOVAR Transcript Curation](TranscriptSplittingAnnovar.md)
- [SnpEff Transcript Curation](TranscriptSplittingSnpeff.md)
- [VEP Transcript Curation](TranscriptSplittingVep.md)
