# ClinVar Data Collection and Curation

## Download ClinVar VCF, variant_summary and variation_allele

Using version 20240107  
Note: variation_allele is provided without archive.
```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2024/clinvar_20240107.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2024/clinvar_20240107.vcf.gz.tbi

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2024-01.txt.gz
gzip -d variant_summary_2024-01.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variation_allele.txt.gz
gzip -d variation_allele.txt.gz
```

## VCF Normalization
Use `bcftools` v1.18

Check the number of variants in the VCF file:

```bash
version=20240107
bcftools stats clinvar_${version}_GRCh38.vcf.gz
```

```
number of records:      2349000
number of no-ALTs:      990
number of SNPs: 2147692
number of MNPs: 7847
number of indels:       186883
number of others:       5588
```

Remove non-standard chromosomes (keep only 1–22, X, Y, MT).  
Reference genome: Ensembl GRCh38 release-111.

```bash
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT  -O z -o clinvar_GRCh38.exChr.vcf.gz clinvar_${version}_GRCh38.vcf.gz

bcftools norm \
  -m -any \
  -O z \
  -cs \
  -f Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
  -o clinvar_GRCh38.norm1.vcf.gz \
  clinvar_GRCh38.exChr.vcf.gz

```

Remove duplicate variants:
```bash
bcftools norm \
  --no-version \
  -d none \
  -O z \
  -o clinvar_GRCh38.norm2.vcf.gz \
  clinvar_GRCh38.norm1.vcf.gz
```

Remove variants with no ALT, non-standard alleles, or length > 50bp:
```bash
bcftools view \
  --no-version \
  --type snps,indels,mnps \
  -e 'ILEN>51 | ILEN<-51 | ALT~"R" | ALT~"Y" | ALT~"M" | ALT~"K" | ALT~"S" | ALT~"W" | ALT~"H" | ALT~"B" | ALT~"V" | ALT~"D" | ALT~"N" | ALT~"*"' \
  -O z -o clinvar_GRCh38.norm3.vcf.gz \
  clinvar_GRCh38.norm2.vcf.gz
```

Count and index the normalized VCF:
```bash
bcftools index clinvar_GRCh38.norm3.vcf.gz
bcftools stats clinvar_GRCh38.norm3.vcf.gz
```
```
number of records:      2338611
number of no-ALTs:      0
number of SNPs: 2147643
number of MNPs: 7846
number of indels:       183122
number of others:       0
```

## Transform VCF to TSV for R or Python
Use `bcftools query`.  
Fix HTML encoding with `sed`:
```bash
bcftools query -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%AF_ESP\t%AF_EXAC\t%AF_TGP\t%ALLELEID\t%CLNDN\t%CLNDNINCL\t%CLNDISDB\t%CLNDISDBINCL\t%CLNHGVS\t%CLNREVSTAT\t%CLNSIG\t%CLNSIGCONF\t%CLNSIGINCL\t%CLNVC\t%CLNVCSO\t%CLNVI\t%DBVARID\t%GENEINFO\t%MC\t%ORIGIN\t%RS\n' clinvar_GRCh38.norm3.vcf.gz > clinvar_GRCh38.norm.tsv

sed -e 's/%3D/=/g' -i clinvar_GRCh38.norm.tsv
```

## variant_summary Curation and Merging
Use R language. See [ClinVar Sample Set Curation](ClinVarSamplesetCuration.md)

## ClinVar CLNREVSTAT 2-star VCF for Annotators
Use R language. See [Create ClinVar VCF](CreateClinVarVCF.md).  
You can use `ClinVar_GRCh38_PLP_Tx1_20240107.vcf` and `ClinVar_GRCh38_BLB_Tx1_20240107.vcf` for the following annotator analysis, or download [ClinVar_GRCh38_PLP_Tx1_20240107.vcf.gz](ClinVar_GRCh38_PLP_Tx1_20240107.vcf.gz) and [ClinVar_GRCh38_BLB_Tx1_20240107.vcf.gz](ClinVar_GRCh38_BLB_Tx1_20240107.vcf.gz) and decompress them.
