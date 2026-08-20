# Rendering and Sharing {#render-share}

This chapter explains how to render and share the tutorial while keeping the computationally heavy methylation analyses separate from the documentation build.

## Render without running heavy code

This is the default.


``` r
bookdown::render_book(
  "index.Rmd",
  output_format = "bookdown::bs4_book",
  params = list(run_code = FALSE)
)
```

This mode is useful for reading, editing, and sharing the documentation. Because the analysis chunks are not executed, readers do not need the IDAT files or the large processed `.rds` objects simply to build the book.

## Running the analysis code

The intended use is to render the book with code execution disabled and run the relevant analysis chunks or companion scripts in the order described by each chapter. Some chunks are teaching templates rather than a complete first-run pipeline.

For a fully automated reproducible analysis, use ordered R scripts or a workflow manager such as `targets`, and keep the book as the explanation and reporting layer. This also avoids hidden dependencies between exploratory script sections.

## Display modes

The HTML book includes a small **Display** control in the lower-right corner. It has two independent settings:

- **Text** switches the page, sidebar, tables, and navigation between light and dark mode.
- **Code** switches inline and fenced code blocks between light and dark mode.

The choices are saved in the browser with `localStorage`, so the same combination is restored when the book is reopened. The small button in the control collapses the panel during reading.

## Recommended sharing bundle

A useful sharing bundle includes:

- `scripts/bookdown_episignatures_tutorial/`;
- the helper functions required by the tutorial;
- a README explaining where public or controlled-access data can be obtained;
- small example phenotype tables that contain no restricted information;
- final CpG lists and trained models when their redistribution is permitted.

Public GEO IDATs can be linked to their repository records rather than duplicated unnecessarily. Controlled or identifiable genomic/clinical data should only be distributed under the permissions that apply to the original study.

## Suggested public README text

```text
This repository contains a tutorial for episignature discovery and validation.
The bookdown source is in scripts/bookdown_episignatures_tutorial/.

To render the documentation without executing analysis code:

  cd scripts/bookdown_episignatures_tutorial
  Rscript render_book.R

To execute the workflow, run the relevant chapter chunks or companion scripts
after installing the required Bioconductor and CRAN packages and placing the
IDAT files in the expected dataset folders.
```

## Session information

For manuscript analyses and public releases, save the software environment used to generate the results:


``` r
sessionInfo()
```

Also record:

- operating system;
- R version;
- Bioconductor version;
- package versions;
- random seeds;
- date of GEO downloads;
- exact CpG list used for each model.

## Final workflow summary

For a new episignature-discovery project:

1. Harmonize and verify phenotype data.
2. Read, normalize, and filter the methylation arrays.
3. Review sample QC, sex checks, probe QC, and technical PCA.
4. Estimate and adjust for blood-cell composition when appropriate.
5. Define discovery and validation samples before feature selection.
6. Use limma or another justified method to identify DMPs in discovery samples only.
7. Select the CpG panel without using validation labels.
8. Train the classifier with fixed preprocessing and predictors.
9. Validate on held-out and, preferably, independent cohorts.
10. Test cross-disease specificity when appropriate.
11. Revalidate before transferring a model to a different array generation or sequencing platform.
12. Record every threshold, random seed, and manual decision.

The two worked examples cover complementary situations: Williams demonstrates de novo signature discovery and independent testing, while Kabuki demonstrates validation of an established signature across cohorts and methylation technologies.

