# HCC_microbiota2022

## Dependencies

- ***R*** (version >= 4.1.2)
- R-packages/extensions
  - dada2 (v. 1.24.0)
  - Biostrings (v. 2.64.0)
  - data.table (v. 1.14.2)
  - vegan (v. 2.6-4)
  - ggplot2 (v. 3.3.6)
  - bit64 (v. 4.0.5)
  - ggtext (v. 0.1.1)
  - egg (v. 0.4.5)
  - ggpubr (v. 0.4.0)
  - DESeq2 (v. 1.36.0)
  - ggrepel (v. 0.9.1)
  - stringr (v. 1.4.0)
  - ComplexHeatmap (v. 2.12.0)
  - circlize (v. 0.4.15)
  - glue (v. 1.6.2)
  - VennDiagram (v. 1.7.3)
  - svglite (v. 2.1.0)
  - topGO (v. 2.48.0)

## Repository structure

- `data/` contains meta info, transcriptome and microbiome data
  - `dada/` results from dada2 processing of 16S data
  - `fastq/` directory for 16S sequencing data in fastq format. Download the fastq files from the ENA-Project with the accession number XYZ and place the files in a subdirectory names "DNA".
  - `meta/` sample information tables
  -  `rnaseq/` transcriptome data (i.e. gene counts)
  - `silva_refs` is a direcotory that ic created as part of the script workflow below. It will contain the 16S taxonomy reference database from SILVA.
- `output/` Directory for data analysis workflow output files
  - `plots`
  - `files` e.g. intermediate files from workflow
- `scripts/` All R-scripts for data analysis.

## Data preprocessing

##### 16S rRNA amplicon data from DNA

The merged fastq files are stored in `data/fastq/DNA/`.

The fastq files are analysed using 'dada' and 'vegan' in the preprocessing script:

```shell
Rscript scripts/00_preprocess_16SDNA.R
```

## Data analysis

##### Cohort overview

A graphical representations of the cohort structure.

```r
# Cohort overview as dot-plot
source("scripts/init_cohort_overview.R")
```

##### Diversity metrics

Alpha diversity

```r
# DNA. Diversity on raw counts
source("scripts/stat_alphaDiversity_16SDNA.R")
```

Beta diversity and NMDS plot

```r
# Beta diversity as Bray-Curtis distance of ASV counts
# Statistics: PERMANOVA as of vegan's adonis()
source("scripts/stat_betaDiversity.R")
```

##### Compositional analysis

Faecal bacteria in blood and liver samples. Definition of <u>faecal bacteria</u>: "The ASV has an relative abundance of >= 1% in at least 1 faecal sample."

```r
# Faecal bacteria proportion in blood and liver
source("scripts/stat_fecalBacteriaProportion_16SDNA.R")
```

Differential analysis of Genus counts using DESeq2

```R
# DNA samples
source("scripts/stat_16S_deseq2_Genus.R")
```

##### Transcriptome (RNA-seq) analysis

Differential gene expression analysis using 'DESeq2'. GO-Term enrichment analysis using 'topGO'

```r
# Main trascriptome analysis
source("scripts/stat_RNA-Seq_deseq2_GOenr2.R")
```

##### RNA-Seq -X- 16S-Liver data co-analysis

*Basic questions: Are there association between the relative abundance of bacterial genera in the liver samples and the transcriptomic profiles?*

deseq analysis design: `design <- ~ genus_abun + age_sc + sex`

where the Genus abundance is normalised using variant-stabilizing transformation as implemented in DESeq2 and susequently z-transformed per sample to speed up convergence in deseq model fitting.

Analysis is limited to HCC samples (n=17).

```r
source("scripts/stat_DNA_Genus_rnaseq_deseq.R")
```

Gene names are colored based on their association with the Condition (see MA-Plots). Red: Up-regulated in HCC, blue: up-regulated in NAFLD.

