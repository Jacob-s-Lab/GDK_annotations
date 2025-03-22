# Create ClinVar P+LP and B+LB VCF
Use R tidyverse Library
```r
library(tidyverse)

load("df_clinvar_GRCh38_BLB_Tx1.RData")
load("df_clinvar_GRCh38_PLP_Tx1.RData")
```

```r
df_clinvar_GRCh38_PLP_Tx1 %>%
  select(chr, start, ref, alt, VariationID, AlleleID) %>%
  # mutate(key = str_c(VariationID, AlleleID, sep = "|")) %>%
  rename(`#CHROM` = chr, POS = start, REF = ref, ALT = alt) %>%
  mutate(ID = ".", QUAL = ".", FILTER = ".") %>%
  mutate(FORMAT = 'GT:AD:DP') %>%
  mutate(ClinVarg37 = "1/1:0,100:100") %>%
  mutate(INFO = str_c("TAG=", key)) %>%
  select(
    `#CHROM`,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    FILTER,
    INFO,
    FORMAT,
    ClinVar
  ) %>%
  write_tsv(file = "ClinVar_GRCh38_PLP_Tx1_20240107.mvcf")

df_clinvar_GRCh38_BLB_Tx1 %>%
  select(chr, start, ref, alt, VariationID, AlleleID) %>%
  # mutate(key = str_c(VariationID, AlleleID, sep = "|")) %>%
  rename(`#CHROM` = chr, POS = start, REF = ref, ALT = alt) %>%
  mutate(ID = ".", QUAL = ".", FILTER = ".") %>%
  mutate(FORMAT = 'GT:AD:DP') %>%
  mutate(ClinVarg37 = "1/1:0,100:100") %>%
  mutate(INFO = str_c("TAG=", key)) %>%
  select(
    `#CHROM`,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    FILTER,
    INFO,
    FORMAT,
    ClinVar
  ) %>%
  write_tsv(file = "ClinVar_GRCh38_BLB_Tx1_20240107.mvcf")
```
Use [vcf header](vcf_header.txt) combind `.mvcf` to `.vcf`
```bash
cat vcf_header.txt ClinVar_GRCh38_PLP_Tx1_20240107.mvcf > ClinVar_GRCh38_PLP_Tx1_20240107.vcf
rm ClinVar_GRCh38_PLP_Tx1_20240107.mvcf

cat vcf_header.txt ClinVar_GRCh38_BLB_Tx1_20240107.mvcf > ClinVar_GRCh38_BLB_Tx1_20240107.vcf
rm ClinVar_GRCh38_BLB_Tx1_20240107.mvcf
```

## Bgzip and Indxing
Use `bcftools` v1.18
```bash
bgzip -k ClinVar_GRCh38_PLP_Tx1_20240107.vcf
bgzip -k ClinVar_GRCh38_BLB_Tx1_20240107.vcf

bcftools index ClinVar_GRCh38_PLP_Tx1_20240107.vcf.gz
bcftools index ClinVar_GRCh38_BLB_Tx1_20240107.vcf.gz
```





