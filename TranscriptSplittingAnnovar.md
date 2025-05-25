After running ANNOVAR annotation, several output files are generated, including `_multianno.txt`, `exonic_variant_function`, and `refGene.variant_function`.

While most users work with the consolidated `_multianno.txt`, here we will process the `exonic_variant_function` and `refGene.variant_function` files separately.

Use the R **tidyverse** library

```r
library(tidyverse)
```

## Processing `refGene.variant_function`

Use the `ClinVar_GRCh38_PLP.refGene.variant_function` or `ClinVar_GRCh38_BLB.refGene.variant_function` file produced by ANNOVAR annotation.

During the ClinVar dataset curation process, the `key` created from `VariationID` and `AlleleID` is also output by ANNOVAR and will serve as a unique identifier.

```r
tsv_all_variant <- read_tsv(
  "your/annovar/output/refGene.variant_function",
  col_names = c("Func_refGene", "Gene_refGene", "no_chr", "no_start", "no_end", "no_ref", "no_alt", "chr", "start", "id", "ref", "alt", "QUAL", "FILTER", "key", "INFO", "FORMAT")
) %>%
  select(key, chr, start, ref, alt, Func_refGene, Gene_refGene)

tsv_all_variant_sepTx <- tsv_all_variant %>%
  mutate(Gene_refGene = str_replace_all(Gene_refGene, '([^\)\d]),', '\1&')) %>%
  separate_rows(Gene_refGene, sep = ",")

rm(tsv_all_variant)
```

## Splitting into intergenic, UTR, splicing, ncRNA_splicing, upstream/downstream, intronic/ncRNA, and exonic categories

```r
tsv_all_variant_intergenic <- tsv_all_variant_sepTx %>%
  filter(Func_refGene == "intergenic") %>%
  rename(Consequence = Func_refGene) %>%
  mutate(NM = str_extract(Gene_refGene, '[^\(]*')) %>%
  mutate(HGVSc = str_extract(Gene_refGene, 'c\.[^\)]*')) %>%
  select(key, chr, start, ref, alt, NM, Consequence, HGVSc)

# UTR regions

tsv_all_variant_utr <- tsv_all_variant_sepTx %>%
  filter(Func_refGene %in% c("UTR5", "UTR3")) %>%
  rename(Consequence = Func_refGene) %>%
  mutate(NM = str_extract(Gene_refGene, '[^\(]*')) %>%
  mutate(HGVSc = str_extract(Gene_refGene, 'c\.[^\)]*')) %>%
  select(key, chr, start, ref, alt, NM, Consequence, HGVSc)

# Splicing categories

tsv_all_variant_sp <- tsv_all_variant_sepTx %>%
  filter(Func_refGene %in% c("splicing", "ncRNA_splicing")) %>%
  rename(Consequence = Func_refGene) %>%
  mutate(NM = str_extract(Gene_refGene, '[^\(]*')) %>%
  mutate(HGVSc = str_extract(Gene_refGene, 'c\.[^\)&]*')) %>%
  select(key, chr, start, ref, alt, NM, Consequence, HGVSc)

# Upstream and downstream regions

tsv_all_variant_stream <- tsv_all_variant_sepTx %>%
  filter(Func_refGene %in% c("upstream", "downstream")) %>%
  rename(Consequence = Func_refGene) %>%
  mutate(NM = str_extract(Gene_refGene, '[^\(]*')) %>%
  mutate(HGVSc = str_extract(Gene_refGene, 'dist=\d+')) %>%
  mutate(HGVSc = if_else(is.na(HGVSc), ".", HGVSc)) %>%
  select(key, chr, start, ref, alt, NM, Consequence, HGVSc)

# Intronic and ncRNA regions

tsv_all_variant_innc <- tsv_all_variant_sepTx %>%
  filter(Func_refGene %in% c("intronic", "ncRNA_intronic", "ncRNA_exonic")) %>%
  rename(Consequence = Func_refGene) %>%
  rename(NM = Gene_refGene) %>%
  select(key, chr, start, ref, alt, NM, Consequence)

# Exonic variants

tsv_all_variant_exonic <- tsv_all_variant_sepTx %>%
  filter(Func_refGene == "exonic")

rm(tsv_all_variant_sepTx)
```

## Processing `refGene.exonic_variant_function`

Use the `ClinVar_GRCh38_PLP.refGene.exonic_variant_function` or `ClinVar_GRCh38_BLB.refGene.exonic_variant_function` file produced by ANNOVAR annotation.

```r
tsv_exon_variant <- read_tsv(
  "your/annovar/output/refGene.exonic_variant_function",
  col_names = c("line", "ExonicFunc_refGene", "AAChange_refGene", "no_chr", "no_start", "no_end", "no_ref", "no_alt", "chr", "start", "id", "ref", "alt", "QUAL", "FILTER", "key", "INFO", "FORMAT")
) %>%
  select(key, chr, start, ref, alt, ExonicFunc_refGene, AAChange_refGene)

tsv_exon_variant_sepTX <- tsv_exon_variant %>%
  separate_rows(AAChange_refGene, sep = ",") %>%
  filter(AAChange_refGene != "")

tsv_exon_variant_sepTXsepHGVS <- tsv_exon_variant_sepTX %>%
  filter(ExonicFunc_refGene != "unknown") %>%
  separate(AAChange_refGene, sep = ":", into = c("Symbol", "Gene_refGene", "exon_NO", "HGVSc", "HGVSp"))

rm(tsv_exon_variant, tsv_exon_variant_sepTX)
```

## Joining exonic variant details

```r
Joint_all_variant_exonic <- tsv_all_variant_exonic %>%
  left_join(
    tsv_exon_variant_sepTXsepHGVS,
    by = c("key", "chr", "start", "ref", "alt", "Gene_refGene")
  )

Joint_all_variant_exonic_unknown <- Joint_all_variant_exonic %>%
  rename(NM = Gene_refGene) %>%
  mutate(
    ExonicFunc_refGene = if_else(is.na(ExonicFunc_refGene), "unknown", ExonicFunc_refGene),
    HGVSc = if_else(is.na(HGVSc), "UNKNOWN", HGVSc)
  ) %>%
  rename(Consequence = ExonicFunc_refGene) %>%
  select(key, chr, start, ref, alt, NM, Consequence, HGVSc, Symbol, HGVSp)

rm(tsv_all_variant_exonic, Joint_all_variant_exonic, tsv_exon_variant_sepTXsepHGVS)
```

## Merging all annotations

```r
df_ClinVar_GRCh38_annovar_refseq <- bind_rows(
  Joint_all_variant_exonic_unknown,
  tsv_all_variant_intergenic,
  tsv_all_variant_innc,
  tsv_all_variant_sp,
  tsv_all_variant_stream,
  tsv_all_variant_utr
) %>%
  mutate(key = str_remove(key, "TAG=")) %>%
  arrange(match(chr, c(1:22, "X", "Y", "MT")), start, ref, alt, NM)

df_ClinVar_GRCh38_annovar_refseq[df_ClinVar_GRCh38_annovar_refseq == ""] <- "."
df_ClinVar_GRCh38_annovar_refseq[is.na(df_ClinVar_GRCh38_annovar_refseq)] <- "."
```

> Note: The processed `df_ClinVar_GRCh38_annovar_refseq` can be used to compare against the `Curated ClinVar TSV` files using the join `key` we defined. The same approach applies to other annotator post-processing outputs.




