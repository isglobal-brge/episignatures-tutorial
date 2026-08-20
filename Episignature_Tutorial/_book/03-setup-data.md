# Setup and Data Organization {#setup-data}

The repository uses a simple organization: one folder per syndrome, one subfolder per dataset, and predictable filenames for phenotype tables, processed methylation objects, figures, and exported results. Keeping this structure consistent makes the later code easier to follow.

## Recommended folder structure

For a discovery dataset:

```text
Williams/GSE119778/
  IDATs/
  Tables/
  Figures/
  tmp_data/
  pheno_GSE119778.csv
  GRset_GSE119778.rds
```

For an external validation dataset:

```text
Williams/GSE66552/
  IDATs/
  Tables/
  Figures/
  pheno_GSE66552.csv
  GRset_GSE66552.rds
  GRset_epi_GSE66552.rds
```

Kabuki follows the same pattern. The long-read and matched-array workbook from the Hildonen et al. study is stored under `Kabuki/PMID41225292/` [@hildonen2026longread].

## Install packages

The workflow uses both Bioconductor and CRAN packages. Installation is normally done once before the analysis, not every time the book is rendered.

The main package groups are:

- **Data access and methylation preprocessing:** `GEOquery`, `minfi`, `FlowSorted.Blood.EPIC`, `maxprobes`, `Biobase`, `S4Vectors`, and `SummarizedExperiment` [@davis2007geoquery; @aryee2014minfi].
- **Differential methylation:** `limma` [@ritchie2015limma].
- **Data handling and plotting:** the `tidyverse` plus the plotting and workbook packages listed in the installation chunk [@wickham2019tidyverse].
- **Machine learning:** `caret` coordinates model training, while packages such as `e1071`, `glmnet`, `ranger`, `kknn`, `gbm`, `xgboost`, and `nnet` provide individual model backends [@kuhn2008caret].
- **Book rendering:** `bookdown` [@xie2016bookdown].

The complete installation list is kept in the code below so readers can install the same dependencies without having to identify them one by one.


``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GEOquery",
  "minfi",
  "limma",
  "FlowSorted.Blood.EPIC",
  "ComplexHeatmap",
  "S4Vectors",
  "SummarizedExperiment",
  "Biobase",
  "maxprobes"
))

install.packages(c(
  "tidyverse",
  "caret",
  "doParallel",
  "e1071",
  "glmnet",
  "ranger",
  "kknn",
  "pls",
  "klaR",
  "gbm",
  "xgboost",
  "nnet",
  "data.table",
  "ggrepel",
  "ggpubr",
  "ggbeeswarm",
  "pals",
  "readxl",
  "openxlsx",
  "RColorBrewer",
  "patchwork",
  "bookdown"
))
```

The analyses in this tutorial use a set of reusable wrapper functions for methylation preprocessing, quality control, machine-learning model training, prediction, and visualization. These functions are loaded from `R/functions_tutorial_epis.R`. Their complete definitions and additional explanations are provided in Chapter \@ref(wrapper-functions).

## Load packages for an analysis session


``` r
library(GEOquery)
library(minfi)
library(tidyverse)
library(limma)
library(caret)
library(doParallel)
library(ggrepel)
library(ggpubr)
library(readxl)
library(openxlsx)

source(file.path("R", "functions_tutorial_epis.R"))

if (file.exists(epiflare_utils)) {
  source(epiflare_utils)
}
```

## Define paths

The examples build paths with `file.path()` instead of hard-coded server locations. This makes it clearer where each input and output lives and makes the workflow easier to move to another computer.


``` r
geo_ID <- "GSE119778"
NDD_syndrome <- "Williams"

path_save <- file.path(project_root, NDD_syndrome, geo_ID)
path_tmp <- file.path(path_save, "tmp_data")
path_idat <- file.path(path_save, "IDATs")
path_figures <- file.path(path_save, "Figures")

ensure_project_dirs(path_save)
```

## Download GEO metadata

GEO datasets do not use one universal set of phenotype column names. The safest approach is to inspect the metadata first and then map the dataset-specific fields into the labels used by this tutorial.


``` r
gse <- GEOquery::getGEO(GEO = geo_ID, destdir = path_idat)[[1]]
geo_pheno_raw <- Biobase::pData(gse)

colnames(geo_pheno_raw)
dplyr::glimpse(geo_pheno_raw)
```

## Harmonized phenotype schema

For combined analyses, the tutorial aims for the following common phenotype schema:

| Column | Meaning |
|--------|---------|
| `title` | GEO sample title or local sample name |
| `geo_accession` | sample key, usually a GSM ID |
| `Sex` | reported or inferred sex as `Male`/`Female` |
| `Case_Cont` | binary label: `Control` or `Case` |
| `Diagnosis` | more detailed label, such as Williams, Dup7, Kabuki, or CHARGE |
| `Age` | age when available; otherwise `NA` |
| `Dataset` | cohort or platform label |

Not every raw parser creates all of these fields immediately. For example, age is unavailable for the Williams discovery cohort, and some dataset labels are added during preprocessing or harmonization. The important point is that the columns required by a later analysis must exist before that analysis is run.

`Case_Cont` is the binary outcome used in the Williams limma model and in classifier training. `Diagnosis` and `Dataset` retain more detailed information for interpretation and plotting.

## Williams phenotype parser

The Williams discovery dataset `GSE119778` comes from the Kimura et al. study [@kimura2020williams]. The script maps its metadata into the tutorial labels as follows:


``` r
pheno_williams <- Biobase::pData(gse) %>%
  dplyr::mutate(
    Sex = dplyr::case_when(
      `gender:ch1` == "Male" ~ "Male",
      `gender:ch1` == "Female" ~ "Female",
      TRUE ~ NA_character_
    ),
    Diagnosis = `diagnosis:ch1`,
    Case_Cont = dplyr::case_when(
      Diagnosis == "Control" ~ "Control",
      TRUE ~ "Case"
    )
  ) %>%
  dplyr::select(title, geo_accession, Sex, Case_Cont, Diagnosis)

readr::write_csv(
  pheno_williams,
  file.path(path_save, paste0("pheno_", geo_ID, ".csv"))
)
```

For the independent Williams validation/comparator dataset `GSE66552`, the diagnosis field is `group:ch1`, and typical controls are labelled `TD control` [@strong2015williams].


``` r
geo_ID_external <- "GSE66552"
path_save_external <- file.path(project_root, "Williams", geo_ID_external)
path_idat_external <- file.path(path_save_external, "IDATs")
ensure_project_dirs(path_save_external)

gse_external <- GEOquery::getGEO(GEO = geo_ID_external, destdir = path_idat_external)[[1]]

pheno_external <- Biobase::pData(gse_external) %>%
  dplyr::mutate(
    Sex = NA_character_,
    Diagnosis = `group:ch1`,
    Case_Cont = dplyr::case_when(
      Diagnosis == "TD control" ~ "Control",
      TRUE ~ "Case"
    )
  ) %>%
  dplyr::select(title, geo_accession, Sex, Case_Cont, Diagnosis)
```

## Kabuki phenotype parsers

The Kabuki workflow combines several datasets, so each one needs its own metadata parser.

For `GSE97362`, which was described by Butcher et al. [@butcher2017chargekabuki]:


``` r
geo_ID <- "GSE97362"
path_save <- file.path(project_root, "Kabuki", geo_ID)
path_idat <- file.path(path_save, "IDATs")
ensure_project_dirs(path_save)

gse <- GEOquery::getGEO(GEO = geo_ID, destdir = path_idat)[[1]]

pheno_gse97362 <- Biobase::pData(gse) %>%
  dplyr::mutate(
    Age = `age (years):ch1`,
    Sex = dplyr::case_when(
      `gender:ch1` == "male" ~ "Male",
      `gender:ch1` == "female" ~ "Female",
      TRUE ~ NA_character_
    ),
    Diagnosis = `disease state:ch1`,
    SampleType = `sample type:ch1`,
    Case_Cont = NA_character_
  ) %>%
  dplyr::select(title, geo_accession, Sex, Case_Cont, Age, Diagnosis, SampleType)
```

For `GSE116300`, described by Sobreira et al. [@sobreira2017kabuki], the script treats both unrelated controls and case parents as `Control` for the binary `Case_Cont` field:


``` r
pheno_gse116300 <- Biobase::pData(gse) %>%
  dplyr::mutate(
    Sex = dplyr::case_when(
      `Sex:ch1` == "Male" ~ "Male",
      `Sex:ch1` == "Female" ~ "Female",
      TRUE ~ NA_character_
    ),
    Diagnosis = `case status:ch1`,
    Case_Cont = dplyr::case_when(
      Diagnosis %in% c("control", "case parent") ~ "Control",
      TRUE ~ "Case"
    ),
    variant_classification = `variant classification:ch1`,
    mutation = `mutation:ch1`
  ) %>%
  dplyr::select(title, geo_accession, Sex, Case_Cont, variant_classification, mutation)
```

Phenotype parsing is not just bookkeeping. An incorrect case/control or diagnosis label can propagate into DMP discovery, model training, and performance assessment. Always inspect the resulting labels with simple tables before continuing.

## Download raw IDATs

Many GEO methylation studies provide raw IDAT files as supplementary archives. The example below shows the download pattern used for `GSE119778`:

```bash
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE119778&format=file" -O GSE119778_RAW.tar
tar -xf GSE119778_RAW.tar -C IDATs
```

After extraction, each array sample should have a matching `*_Grn.idat` and `*_Red.idat` file. `minfi::read.metharray.exp()` reads these paired red and green channel files into the raw methylation object [@aryee2014minfi].



