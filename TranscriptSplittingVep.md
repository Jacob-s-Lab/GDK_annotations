# VCF to TSV

After running the VEP annotation, we will obtain a VCF file corresponding to our input. To facilitate downstream processing, we convert it to a TSV format using `bcftools`. During this conversion, the TAG used as the `key` is extracted along with the `CSQ` field.

Use `bcftools` v1.18 and `split-vep` plugin:
```bash
bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%TAG\t%CSQ\n' -A tab -d  /your/vep/output/vep.vcf > /your/vep/output/vep.tsv
```

# Processing `vep.tsv`

We’ll use the R **tidyverse** package for processing:

```r
library(tidyverse)

tsv_query <- read_tsv(
  file("/your/vep/output/vep.tsv"), 
  col_names = c(
    "chr",
    "start",
    "ref",
    "alt",
    "key",
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "VARIANT_CLASS",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "CANONICAL",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "NP",
    "REFSEQ_MATCH",
    "REFSEQ_OFFSET",
    "GIVEN_REF",
    "USED_REF",
    "BAM_EDIT",
    "HGVS_OFFSET",
    "HGVSg"
  ),
  col_types = cols(start = col_number(), .default = col_character())
)

tsv_nm <- tsv_query %>%
  select(
    key, chr, start, ref, alt,
    VARIANT_CLASS, SYMBOL, Consequence, IMPACT,
    Gene, Feature, Feature_type, BIOTYPE, STRAND, FLAGS,
    EXON, INTRON, cDNA_position, CDS_position, DISTANCE,
    NP, Protein_position, Amino_acids, Codons,
    HGVSg, HGVSc, HGVSp,
    CANONICAL, MANE_SELECT, MANE_PLUS_CLINICAL, HGVS_OFFSET
  ) %>%
  mutate(NM = str_extract(Feature, '[^\\.]+'))

rm(tsv_query)
```

## Retain Valid RefSeq Transcripts

```r
df_ClinVar_GRCh38_vep_refseq <- df_ClinVar_GRCh38_PLP_vep_refseq %>%
  filter(str_detect(NM, '^N\\w_[0-9]+'))

rm(tsv_nm)
```

> **Note:** The processed `df_ClinVar_GRCh38_vep_refseq` can be used to compare against the curated ClinVar TSV files by joining on the `key` field. The same approach applies to other annotator post-processing outputs.
