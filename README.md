# Compare Three Annotators For Genetic Interpretation

---

# ClinVar Data Collection and Curation

## Download ClinVar VCF & variant_summary

Use Version 20240107

```bash
mkdir ClinVar_v20240107
cd ClinVar_v20240107

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2024/clinvar_20240107.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2024/clinvar_20240107.vcf.gz.tbi

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/gene_specific_summary_2024-01.txt.gz
```

## VCF Normalization
Use `bcftools` v1.18

Check VCF Variant Count

```bash
cd ClinVar_v20240107

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

Remove non 1~22, X, Y, MT chrmosome
Using Ensembl GRCh38 release-111 Reference Genome

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

Remove duplicated Variant
```bash
bcftools norm \
  --no-version \
  -d none \
  -O z \
  -o clinvar_GRCh38.norm2.vcf.gz \
  clinvar_GRCh38.norm1.vcf.gz
```

Remove no-ALTs, others and > 50bp Variants
```bash
bcftools view \
  --no-version \
  --type snps,indels,mnps \
  -e 'ILEN>51 | ILEN<-51 | ALT~"R" | ALT~"Y" | ALT~"M" | ALT~"K" | ALT~"S" | ALT~"W" | ALT~"H" | ALT~"B" | ALT~"V" | ALT~"D" | ALT~"N" | ALT~"*"' \
  -O z -o clinvar_GRCh38.norm3.vcf.gz \
  clinvar_GRCh38.norm2.vcf.gz
```

Count and Index Normalized VCF
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

## Transform vcf to tsv file for R or Python
Use `bcftools` Query
fix HTML encode with `sed`
```bash
bcftools query -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%AF_ESP\t%AF_EXAC\t%AF_TGP\t%ALLELEID\t%CLNDN\t%CLNDNINCL\t%CLNDISDB\t%CLNDISDBINCL\t%CLNHGVS\t%CLNREVSTAT\t%CLNSIG\t%CLNSIGCONF\t%CLNSIGINCL\t%CLNVC\t%CLNVCSO\t%CLNVI\t%DBVARID\t%GENEINFO\t%MC\t%ORIGIN\t%RS\n' clinvar_GRCh38.norm3.vcf.gz > clinvar_GRCh38.norm.tsv

sed -e 's/%3D/=/g' -i clinvar_GRCh38.norm.tsv
```


