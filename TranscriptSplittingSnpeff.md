# VCF to TSV

After running the SnpEff annotation, we will obtain a VCF file corresponding to our input. To facilitate downstream processing, we convert it to a TSV format using `bcftools`. During this conversion, the TAG used as the `key` is extracted along with the `ANN` field; you may optionally exclude the `LOF` and `NMD` columns as needed.

Use `bcftools` v1.18:
```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TAG\t%ANN\t%LOF\t%NMD\n' /your/snpeff/output/snpeff.vcf > /your/snpeff/output/snpeff.tsv
```

# Processing `snpeff.tsv`

We’ll use the R **tidyverse** package for processing:

```r
library(tidyverse)

tsv_ann <- read_tsv(
  "/your/snpeff/output/snpeff.tsv",
  col_names = c("chr", "start", "ref", "alt", "key", "ANN", "LOF", "NMD")
)

tsv_ann_sep1 <- tsv_ann %>%
  separate_rows(ANN, sep = ",")

tsv_ann_sep2 <- tsv_ann_sep1 %>%
  separate(
    ANN,
    c(
      "Allele",
      "Annotation",
      "Annotation_Impact",
      "Gene_Name",
      "Gene_ID",
      "Feature_Type",
      "Feature_ID",
      "Transcript_BioType",
      "Rank",
      "HGVSc",
      "HGVSp",
      "cDNA",
      "CDS",
      "AA",
      "Distance",
      "INFO"
    ),
    sep = "\\|", remove = TRUE, convert = FALSE
  )

# Replace empty strings with a dot

tsv_ann_sep2[tsv_ann_sep2 == ""] <- "."

tsv_ann_sep3 <- tsv_ann_sep2 %>%
  select(
    key, chr, start, ref, alt,
    Gene_Name, Gene_ID, Annotation, Annotation_Impact,
    Feature_ID, Feature_Type, Transcript_BioType,
    Rank, cDNA, CDS, AA, Distance,
    HGVSc, HGVSp,
    LOF, NMD
  ) %>%
  mutate(NM = str_extract(Feature_ID, '[^\\.]+'))

rm(tsv_ann, tsv_ann_sep1, tsv_ann_sep2)
```

## Retain Valid RefSeq Transcripts

```r
df_ClinVar_GRCh38_snpeff_refseq <- tsv_ann_sep3 %>%
  filter(str_detect(Feature_ID, '^N\\w_[0-9]+\\.[0-9]+'))

rm(tsv_ann_sep3)
```

> **Note:** The processed `df_ClinVar_GRCh38_snpeff_refseq` can be used to compare against the curated ClinVar TSV files by joining on the `key` field. The same approach applies to other annotator post-processing outputs.
