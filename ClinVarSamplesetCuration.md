# Load variation_allele.txt

Use R tidyverse Library
```r
library(tidyverse)
```

Just Keep `Variant Type` Allele
```r
tsv_variation_allele <- read_tsv(
  "variation_allele.txt",
  skip = 8,
  col_names = c("VariationID", "Type", "AlleleID", "Interpreted")
)

tsv_variation_allele_filter <- tsv_variation_allele %>%
  filter(Interpreted == "yes") %>%
  filter(Type == "Variant")

rm(tsv_variation_allele)
```

# Load variant_summary.txt
```r
tsv_variant_summary <- read_tsv("variant_summary_2024-01.txt")
```
4,804,422 Variant

## Filter out non `Variant Type` Allele
```r
tsv_variant_summary_filter <- tsv_variant_summary  %>%
  filter(`#AlleleID` %in% tsv_variation_allele_filter$AlleleID)

rm(tsv_variation_allele_filter, tsv_variant_summary)
```

## Filter out GRCh38 and Curation
```r
vsum_clinvar_g38 <- tsv_variant_summary_filter %>%
  filter(Assembly == "GRCh38") %>%
  filter(!(ReferenceAlleleVCF %in% c("na", "Un"))) %>%
  filter(!(AlternateAlleleVCF %in% c("na", "Un"))) %>%
  filter(!(Chromosome %in% c("na", "Un"))) %>%
  rename(AlleleID = `#AlleleID`)

rm(tsv_variant_summary_filter)
```

---

# Merge VCF and variant_summary

## Load VCF Normalized Tsv
```r
library(tidyverse)

vcf_clinvar_g38 <- read_tsv(
  "clinvar_GRCh38.norm.tsv",
  col_names = c(
    "chr", "start", "end", "ref", "alt",
    "AF_ESP", "AF_EXAC", "AF_TGP", "AlleleID",
    "CLNDN", "CLNDNINCL", "CLNDISDB", "CLNDISDBINCL", "CLNHGVS", "CLNREVSTAT", "CLNSIG", "CLNSIGCONF", "CLNSIGINCL", "CLNVC", "CLNVCSO", "CLNVI",
    "DBVARID", "GENEINFO", "MC", "ORIGIN", "RS"
  ),
  na = ".",
  col_types = cols(
    chr = col_character(),
    ref = col_character(),
    alt = col_character(),
    AF_ESP = col_character(),
    AF_EXAC = col_character(),
    AF_TGP = col_character(),
    CLNDN = col_character(),
    CLNDNINCL = col_character(),
    CLNDISDB = col_character(),
    CLNDISDBINCL = col_character(),
    CLNHGVS = col_character(),
    CLNREVSTAT = col_character(),
    CLNSIG = col_character(),
    CLNSIGCONF = col_character(),
    CLNSIGINCL = col_character(),
    CLNVC = col_character(),
    CLNVCSO = col_character(),
    CLNVI = col_character(),
    DBVARID = col_character(),
    RS  = col_character()
  )
)
```

## VCF CLNSIG Curation
```r
vcf_clinvar_g38 <- vcf_clinvar_g38 %>%
  mutate(
    ClinVar_CLNSIG = case_when(
      (str_detect(CLNSIG, "Conflicting") & str_detect(CLNSIGCONF, "athogenic")) ~ "Conflicting_PLP",
      str_detect(CLNSIG, "Conflicting") ~ "Conflicting",
      str_detect(CLNSIG, "Pathogenic") ~ "Pathogenic",
      str_detect(CLNSIG, "Likely_pathogenic") ~ "Likely_pathogenic",
      str_detect(CLNSIG, "Benign") ~ "Benign",
      str_detect(CLNSIG, "Likely_benign") ~ "Likely_benign",
      str_detect(CLNSIG, "drug_response") ~ "drug_response",
      str_detect(CLNSIG, "risk_factor") ~ "risk_factor",
      str_detect(CLNSIG, "not_provided") ~ "not_provided",
      str_detect(CLNSIG, "Uncertain_significance") ~ "Uncertain_significance",
      (!is.na(CLNSIG)) ~ "other"
    )
  )
```

## Extract variant_summary `VariationID`, `Type`, `Name` ... to VCF
Join by `chr` and `AlleleID`
```r
df_clinvar_GRCh38 <- vcf_clinvar_g38 %>%
  select(1:5, AlleleID, ClinVar_CLNSIG, CLNSIGCONF, CLNREVSTAT, MC) %>%
  left_join(
    vsum_clinvar_g38 %>%
      rename(chr = Chromosome) %>%
      select(chr, AlleleID, VariationID, Type, Name, GeneSymbol, RCVaccession, PhenotypeIDS, PhenotypeList, NumberSubmitters, OriginSimple),
    by = c("chr", "AlleleID")
  )

save(df_clinvar_GRCh38, file = "df_clinvar_GRCh38_20240107.RData")
rm(vcf_clinvar_g37, vcf_clinvar_g38, vsum_clinvar_g37, vsum_clinvar_g38)
```

---

# ClinVar CLNREVSTAT 2 star Variants Curation

## Remove MT and variant without MC or Multi-MC
```r
df_clinvar_GRCh38_PLP_sep <- df_clinvar_GRCh38 %>%
  filter(ClinVar_CLNSIG %in% c("Pathogenic", "Likely_pathogenic")) %>%
  filter(CLNREVSTAT %in% c("criteria_provided,_multiple_submitters,_no_conflicts", "reviewed_by_expert_panel", "practice_guideline")) %>%
  filter(chr != "MT") %>%
  select(chr:ClinVar_CLNSIG, CLNREVSTAT, MC, VariationID, Type, Name, GeneSymbol) %>%
  filter(!is.na(MC)) %>%
  mutate(NM = str_extract(Name, 'N\\w_\\d+')) %>%
  filter(!is.na(NM)) %>%
  filter(!is.na(VariationID)) %>%
  filter(!str_detect(MC, ",")) %>%
  mutate(key = str_c(VariationID, "|", AlleleID)) %>%
  mutate(HGVSc = str_extract(Name, 'c\\.[^\\s\\(]+')) %>%
  mutate(HGVSp = str_extract(Name, '\\((p\\..*)\\)$', group = 1)) %>%
  mutate(Consequence = str_extract(MC, '[^\\|]+$')) %>%
  select(key, chr:ClinVar_CLNSIG, CLNREVSTAT, Consequence, VariationID, Type, Name, NM, GeneSymbol, HGVSc, HGVSp)

df_clinvar_GRCh38_BLB_sep <- df_clinvar_GRCh38 %>%
  filter(ClinVar_CLNSIG %in% c("Benign", "Likely_benign")) %>%
  filter(CLNREVSTAT %in% c("criteria_provided,_multiple_submitters,_no_conflicts", "reviewed_by_expert_panel", "practice_guideline")) %>%
  filter(chr != "MT") %>%
  select(chr:ClinVar_CLNSIG, CLNREVSTAT, MC, VariationID, Type, Name, GeneSymbol) %>%
  filter(!is.na(MC)) %>%
  mutate(NM = str_extract(Name, 'N\\w_\\d+')) %>%
  filter(!is.na(NM)) %>%
  filter(!is.na(VariationID)) %>%
  filter(!str_detect(MC, ",")) %>%
  mutate(key = str_c(VariationID, "|", AlleleID)) %>%
  mutate(HGVSc = str_extract(Name, 'c\\.[^\\s\\(]+')) %>%
  mutate(HGVSp = str_extract(Name, '\\((p\\..*)\\)$', group = 1)) %>%
  mutate(Consequence = str_extract(MC, '[^\\|]+$')) %>%
  select(key, chr:ClinVar_CLNSIG, CLNREVSTAT, Consequence, VariationID, Type, Name, NM, GeneSymbol, HGVSc, HGVSp)

```

## Consequence Normalization

```r
df_clinvar_GRCh38_PLP_Tx1 <- df_clinvar_GRCh38_PLP_sep %>%
  mutate(
    Consequence = case_when(
      
      Consequence == "splice_acceptor_variant" ~ "splicing_variant",
      Consequence == "splice_donor_variant" ~ "splicing_variant",
      
      str_detect(Consequence, "inframe_") ~ "inframe_variant",
      
      Consequence == "non-coding_transcript_variant" ~ "non_coding_transcript_variant",
      
      TRUE ~ Consequence
      
    )
  )

df_clinvar_GRCh38_BLB_Tx1 <- df_clinvar_GRCh38_BLB_sep %>%
  mutate(
    Consequence = case_when(
      
      Consequence == "splice_acceptor_variant" ~ "splicing_variant",
      Consequence == "splice_donor_variant" ~ "splicing_variant",
      
      str_detect(Consequence, "inframe_") ~ "inframe_variant",
      
      Consequence == "non-coding_transcript_variant" ~ "non_coding_transcript_variant",
      
      TRUE ~ Consequence
      
    )
  )

rm(df_clinvar_GRCh38_PLP_sep)
rm(df_clinvar_GRCh38_BLB_sep)

```

## ClinVar HGVSc Curation

```r
df_clinvar_GRCh38_BLB_Tx1 <- df_clinvar_GRCh38_BLB_Tx1 %>%
  mutate(HGVSc = if_else(str_detect(HGVSc, 'del[ATCG]+'), str_replace(HGVSc, '(.*del)[ATCG]+', '\\1'), HGVSc)) %>%
  mutate(HGVSc = if_else(str_detect(HGVSc, 'delins[ATCG]+'), str_replace(HGVSc, '(.*delins)[ATCG]+', '\\1'), HGVSc)) %>%
  mutate(HGVSc = if_else(str_detect(HGVSc, 'dup[ATCG]+'), str_replace(HGVSc, '(.*dup)[ATCG]+', '\\1'), HGVSc)) %>%
  mutate(HGVSc = if_else(str_detect(HGVSc, 'ins[ATCG]+'), str_replace(HGVSc, '(.*ins)[ATCG]+', '\\1'), HGVSc)) %>%
  mutate(HGVSc = if_else(str_detect(HGVSc, 'inv[ATCG]+'), str_replace(HGVSc, '(.*inv)[ATCG]+', '\\1'), HGVSc))

df_clinvar_GRCh38_PLP_Tx1 <- df_clinvar_GRCh38_PLP_Tx1 %>%
  mutate(HGVSc = if_else(str_detect(HGVSc, 'del[ATCG]+'), str_replace(HGVSc, '(.*del)[ATCG]+', '\\1'), HGVSc)) %>%
  mutate(HGVSc = if_else(str_detect(HGVSc, 'delins[ATCG]+'), str_replace(HGVSc, '(.*delins)[ATCG]+', '\\1'), HGVSc)) %>%
  mutate(HGVSc = if_else(str_detect(HGVSc, 'dup[ATCG]+'), str_replace(HGVSc, '(.*dup)[ATCG]+', '\\1'), HGVSc)) %>%
  mutate(HGVSc = if_else(str_detect(HGVSc, 'ins[ATCG]+'), str_replace(HGVSc, '(.*ins)[ATCG]+', '\\1'), HGVSc)) %>%
  mutate(HGVSc = if_else(str_detect(HGVSc, 'inv[ATCG]+'), str_replace(HGVSc, '(.*inv)[ATCG]+', '\\1'), HGVSc))
```

```r
save(df_clinvar_GRCh38_BLB_Tx1, file = "df_clinvar_GRCh38_BLB_Tx1.RData")
save(df_clinvar_GRCh38_PLP_Tx1, file = "df_clinvar_GRCh38_PLP_Tx1.RData")
```




