---
title: "Episignatures Detection: A Practical Tutorial"
subtitle: "From methylation arrays to CpG panels, classifiers, and external validation"
author: "Albert Alegret-Garcia, Alejandro Caceres, Luis A. Perez-Jurado, Juan R. Gonzalez | Bioinformatics Research Group in Genetic Epidemiology, ISGlobal"
date: "2026-08-19"
site: bookdown::bookdown_site
documentclass: book
bibliography: references.bib
link-citations: yes
colorlinks: yes
params:
  run_code: false
  project_root: ".."
  epiflare_utils: "/PROJECTES/GENOMICS/aalegret/Epimutations_Episignatures/EpiFLARE/R/utils.R"
  discovery_geo: "GSE119778"
  discovery_syndrome: "Williams"
  external_williams_geo: "GSE66552"
  kabuki_training_geo: "GSE97362"
  kabuki_long_read_id: "PMID41225292"
output:
  bookdown::bs4_book:
    split_by: chapter
    theme:
      primary: "#0f766e"
    css: assets/book-style.css
    includes:
      in_header: assets/header.html
  bookdown::pdf_book:
    toc: true
    toc_depth: 3
    number_sections: true
---



# Preface {-}

This book is a practical tutorial for discovering and testing DNA methylation episignatures in rare disease. It is written for readers with different levels of computational experience: bioinformaticians can follow the analysis structure, while geneticists and other non-specialist readers can use the text to understand what each step is doing and why it matters.

It was developed by the **Bioinformatics Research Group in Genetic Epidemiology, ISGlobal** as a practical companion to the manuscript chapter on episignatures and machine learning in rare diseases.

The tutorial follows two complementary worked examples:

1. **Williams syndrome — de novo discovery:** `scripts/Williams_analysis.R` uses `GSE119778` to identify DMPs, select a candidate episignature CpG panel, train classifiers, and evaluate samples that were kept out of discovery. The independent `GSE66552` cohort is then used for external testing [@kimura2020williams; @strong2015williams].
2. **Kabuki syndrome — existing signature validation:** `scripts/Kabuki/` starts from an established Kabuki CpG panel, trains the tutorial SVM using confirmed samples from `GSE97362`, evaluates independent array cohorts (`GSE116300`, `GSE218186`), and tests whether the disease-associated signal remains detectable in ONT and PacBio methylation profiles from the long-read study [@butcher2017chargekabuki; @sobreira2017kabuki; @hildonen2026longread].

The code chunks are disabled by default so the book can be rendered without downloading IDAT files or loading large `.rds` objects. The chapters are intended to be read in order and the relevant code run interactively or through the companion scripts. Some scripts originated as exploratory analyses, so the tutorial emphasizes dependencies, input checks, and the order in which discovery and validation steps should occur.

This is a research and teaching workflow, not a clinical diagnostic report. Episignature interpretation requires appropriate reference cohorts, technical QC, independent validation, cross-disease specificity assessment where possible, and careful consideration of platform and biological limitations [@arefeshghi2020episignatures; @kerkhof2024reporting].

## How to render

From R:


``` r
setwd("scripts/bookdown_episignatures_tutorial")
bookdown::render_book(
  "index.Rmd",
  output_format = "bookdown::bs4_book",
  params = list(run_code = FALSE)
)
```

Or run the small wrapper:

```bash
cd scripts/bookdown_episignatures_tutorial
Rscript render_book.R
```

## Learning goals

By the end of the tutorial, readers will be able to:

- explain the difference between a localized methylation abnormality and a multivariate disease-associated episignature;
- recognize the main objects produced during Illumina methylation-array preprocessing;
- follow the sequence from sample/probe QC to DMP discovery, CpG selection, model training, and independent testing;
- understand why feature selection must be kept separate from validation data;
- distinguish de novo episignature discovery from testing an established CpG panel;
- interpret classifier results together with PCA, disease comparators, and platform effects.

## Repository paths used in the examples

The book assumes it is rendered from `scripts/bookdown_episignatures_tutorial/`.


``` r
project_root
path_williams
path_kabuki
path_discovery
path_external_williams
path_kabuki_training
path_kabuki_long_read
```

Outside the original server environment, set `params$project_root` to the folder that contains the `Williams/` and `Kabuki/` analysis directories. The tutorial then constructs dataset-specific paths from that root rather than relying on the original absolute HPC paths.



