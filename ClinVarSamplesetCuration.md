# Load Normalized Tsv

Use R tidyverse Library
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

# ClinVar 2 star Variants Curation

Remove MT and variant without MC or Multi-MC
```r
df_clinvar_GRCh38_PLP_sep <- df_clinvar_GRCh38 %>%
  filter(ClinVar_CLNSIG %in% c("Pathogenic", "Likely_pathogenic")) %>%
  filter(CLNREVSTAT %in% c("criteria_provided,_multiple_submitters,_no_conflicts", "reviewed_by_expert_panel", "practice_guideline")) %>%
  # 58,843
  filter(chr != "MT") %>%
  # 58,742 (-101)
  select(chr:ClinVar_CLNSIG, CLNREVSTAT, MC, VariationID, Type, Name, GeneSymbol) %>%
  filter(!is.na(MC)) %>%
  # 58,474 (-269)
  mutate(NM = str_extract(Name, 'N\\w_\\d+')) %>%
  filter(!is.na(NM)) %>%
  # 58,473 (-1)
  filter(!is.na(VariationID)) %>%
  # 58,473 (-0) %>%
  filter(!str_detect(MC, ",")) %>%
  # 44,351 (-14122)
  mutate(key = str_c(VariationID, "|", AlleleID)) %>%
  mutate(HGVSc = str_extract(Name, 'c\\.[^\\s\\(]+')) %>%
  # mutate(Symbol = str_extract(Name, '\\((.*)\\)[^$]', group = 1)) %>%
  mutate(HGVSp = str_extract(Name, '\\((p\\..*)\\)$', group = 1)) %>%
  mutate(Consequence = str_extract(MC, '[^\\|]+$')) %>%
  select(key, chr:ClinVar_CLNSIG, CLNREVSTAT, Consequence, VariationID, Type, Name, NM, GeneSymbol, HGVSc, HGVSp)

df_clinvar_GRCh38_BLB_sep <- df_clinvar_GRCh38 %>%
  filter(ClinVar_CLNSIG %in% c("Benign", "Likely_benign")) %>%
  filter(CLNREVSTAT %in% c("criteria_provided,_multiple_submitters,_no_conflicts", "reviewed_by_expert_panel", "practice_guideline")) %>%
  # 147,544
  filter(chr != "MT") %>%
  # 147,479 (-65)
  select(chr:ClinVar_CLNSIG, CLNREVSTAT, MC, VariationID, Type, Name, GeneSymbol) %>%
  filter(!is.na(MC)) %>%
  # 146,891 (-588)
  mutate(NM = str_extract(Name, 'N\\w_\\d+')) %>%
  filter(!is.na(NM)) %>%
  # 146,875 (-16)
  filter(!is.na(VariationID)) %>%
  # 146,875 (-0) %>%
  filter(!str_detect(MC, ",")) %>%
  # 120,198 (-26677)
  mutate(key = str_c(VariationID, "|", AlleleID)) %>%
  mutate(HGVSc = str_extract(Name, 'c\\.[^\\s\\(]+')) %>%
  # mutate(Symbol = str_extract(Name, '\\((.*)\\)[^$]', group = 1)) %>%
  mutate(HGVSp = str_extract(Name, '\\((p\\..*)\\)$', group = 1)) %>%
  mutate(Consequence = str_extract(MC, '[^\\|]+$')) %>%
  select(key, chr:ClinVar_CLNSIG, CLNREVSTAT, Consequence, VariationID, Type, Name, NM, GeneSymbol, HGVSc, HGVSp)

```

# Consequence Normalization

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

