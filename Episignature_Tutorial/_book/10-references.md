# References and Citation Guide {#references .unnumbered}

This chapter explains how scientific and software references are handled in the tutorial. All citation records are stored in `references.bib`, and the `.Rmd` chapters refer to them by citation key.

## How to cite inside chapters

Pandoc/bookdown citations use square brackets:

```text
DNA methylation arrays were processed with minfi [@aryee2014minfi].
Differential methylation was modelled with limma [@ritchie2015limma].
```

Multiple citations go in the same brackets:

```text
The workflow uses minfi, limma, and caret [@aryee2014minfi; @ritchie2015limma; @kuhn2008caret].
```

Narrative citations use `@key`:

```text
@fortin2014funnorm introduced functional normalization for 450K methylation arrays.
```

## Citation checklist

A complete episignature tutorial should cite the source of both the data and the method. For this project, that includes:

- the original studies associated with the Williams and Kabuki datasets, together with GEO records when useful [@kimura2020williams; @strong2015williams; @butcher2017chargekabuki; @sobreira2017kabuki];
- the publication describing the fixed Kabuki episignature used for validation [@arefeshghi2017kabuki; @arefeshghi2020episignatures];
- the long-read Kabuki study when ONT, PacBio, or matched array data are discussed [@hildonen2026longread];
- methodological references for avoiding selection bias when feature selection and classifier validation are combined [@ambroise2002selection];
- `minfi` for methylation-array processing [@aryee2014minfi];
- functional normalization when used [@fortin2014funnorm];
- `limma` for differential methylation modelling [@ritchie2015limma];
- beta/M-value methodology when those representations are explained [@du2010comparison];
- `caret` and the relevant modelling backends for classifier training [@kuhn2008caret; @e1071; @friedman2010glmnet; @wright2017ranger; @kknn; @mevik2007pls; @klar; @greenwell2022gbm; @chen2016xgboost; @venables2002modern];
- blood cell-mixture estimation methods when cell proportions are modelled [@houseman2012cellmixture; @salas2018blood];
- software used to create figures, tables, or the book itself when those outputs are redistributed [@wickham2019tidyverse; @ggrepel; @ggpubr; @ggbeeswarm; @rcbrewer; @patchwork; @readxl; @openxlsx; @xie2016bookdown].

The bibliography also retains GEO records for `GSE119778`, `GSE66552`, `GSE97362`, `GSE116300`, and `GSE218186` [@geoGSE119778; @geoGSE66552; @geoGSE97362; @geoGSE116300; @geoGSE218186]. In the scientific chapters, the peer-reviewed study is cited when a claim depends on the published cohort or result, while the GEO entry can be used to identify the downloadable dataset.

## Updating package citations

R packages can change their recommended citation as versions are updated. The included `references.bib` contains working citations for the packages used in this tutorial, but these should be checked against the R environment used for the final public release. `knitr::write_bib()` can generate version-specific entries:


``` r
knitr::write_bib(
  c(
    "base",
    "BiocManager",
    "GEOquery",
    "minfi",
    "limma",
    "FlowSorted.Blood.EPIC",
    "ComplexHeatmap",
    "S4Vectors",
    "SummarizedExperiment",
    "Biobase",
    "maxprobes",
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
  ),
  file = "references_packages.bib"
)
```

Before release, compare `references_packages.bib` with the tutorial bibliography and update package entries when the recommended citation has changed. Keep the scientific article references separate from these version-dependent software citations.

## How to add a new reference

Add new scientific or software references to `references.bib` with a unique citation key. For example:

```bibtex
@article{myStudy2026,
  title   = {Title of the study},
  author  = {Surname, Name and Surname, Name},
  journal = {Journal Name},
  year    = {2026},
  volume  = {1},
  number  = {1},
  pages   = {1--10},
  doi     = {10.0000/example}
}
```

The new key can then be cited in any `.Rmd` file:

```text
This cohort was described previously [@myStudy2026].
```

During rendering, bookdown resolves the citation keys and prints the formatted bibliography below. If a citation appears as an unresolved key, first check that the key exists in `references.bib` and is spelled exactly the same way in the `.Rmd` file.

## Reference list

<div id="refs"></div>



