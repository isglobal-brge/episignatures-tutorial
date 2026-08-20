# Preprocessing and Quality Control {#preprocessing-qc}

This chapter explains the main preprocessing and quality-control steps used for the Williams and Kabuki methylation-array datasets. The goal is to show what each check is doing before the data are used for DMP discovery or classification.

## Raw object: `RGChannelSet`

IDATs are first read into an `RGChannelSet`, which stores the raw red and green channel intensities before normalization [@aryee2014minfi]. In a laboratory workflow, a sample sheet can provide sample identifiers and array positions. For these public GEO examples, the script extracts the GEO accession, slide, and array position from the IDAT filenames instead.


``` r
RGset <- minfi::read.metharray.exp(base = path_idat)
RGset
```

The IDAT filename downloaded from GEO starts with the GSM accession, followed by slide and array position. In many local laboratory datasets, filenames contain only the slide and array position:

```text
# GEO
GSM3383551_9379086005_R01C01_Grn.idat

# In-house arrays.
9379086005_R01C01_Grn.idat
```

## Sample detection P-value filter

Detection P-values measure whether probe signal can be distinguished from background. In this repository, the first coarse sample-level filter calculates the **mean** detection P-value across probes for each sample and keeps samples with a mean below 0.05.

This is the rule implemented by these scripts; it is not a universal methylation-array QC threshold. Other workflows may use the fraction of failed probes or additional intensity-based metrics, so the exact rule should be reported when results are published.


``` r
detP <- minfi::detectionP(RGset)
sample_detP <- colMeans(detP, na.rm = TRUE)
print(table(sample_detP < 0.05))

RGset <- RGset[, sample_detP < 0.05]
```
This is only the first sample-level check. Metadata matching, sex prediction, and PCA are reviewed later before samples are used for discovery or training.

## Attach phenotype metadata

Because these examples are read without a sample sheet, the phenotype table is attached after the IDATs are loaded. The GSM accession is extracted from the filename and used as the key that links each array to its metadata.


``` r
pheno <- readr::read_csv(
  file.path(path_save, paste0("pheno_", geo_ID, ".csv")),
  show_col_types = FALSE
)

colData(RGset) <- data.frame(
  geo_accession = sub("_.*", "", colnames(RGset)),
  idat_name = colnames(RGset),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(pheno, by = "geo_accession") %>%
  tibble::column_to_rownames("geo_accession") %>%
  S4Vectors::DataFrame()
```

A simple check confirms that all loaded samples received a case/control label:


``` r
stopifnot(!any(is.na(RGset$Case_Cont)))
table(RGset$Case_Cont, useNA = "ifany")
```

## Sex prediction

Sex prediction provides an additional sample-level QC check. A disagreement between reported and methylation-predicted sex can indicate a sample-label problem, a metadata error, or another issue that deserves review. A discordant result is therefore a flag to investigate, not an automatic reason to discard the sample.


``` r
MSet <- minfi::preprocessNoob(RGset)
MSet <- minfi::mapToGenome(MSet)

predSex <- as.data.frame(minfi::getSex(MSet)) %>%
  tibble::rownames_to_column(var = "geo_accession") %>%
  dplyr::mutate(
    predSex = dplyr::case_when(
      predictedSex == "M" ~ "Male",
      predictedSex == "F" ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(geo_accession, predSex)

colData(RGset) <- as.data.frame(SummarizedExperiment::colData(RGset)) %>%
  tibble::rownames_to_column("geo_accession") %>%
  dplyr::left_join(predSex, by = "geo_accession") %>%
  dplyr::mutate(
    Concordant_Sex = dplyr::case_when(
      is.na(Sex) ~ NA_character_,
      Sex == predSex ~ "Concordant",
      TRUE ~ "Discordant"
    )
  ) %>%
  tibble::column_to_rownames("geo_accession") %>%
  S4Vectors::DataFrame()

table("Reported" = RGset$Sex, "Predicted" = RGset$predSex, useNA = "ifany")
table(RGset$Concordant_Sex, useNA = "ifany")
```

When reported sex is unavailable, concordance cannot be checked. The predicted value can still be stored for plotting or, if scientifically justified, used as a covariate in later analyses.

## Normalization

The Williams discovery workflow uses functional normalization, a method developed for Illumina methylation arrays [@fortin2014funnorm]:


``` r
GRset <- minfi::preprocessFunnorm(RGset, sex = RGset$Sex)
```

The Kabuki datasets are not all processed in exactly the same way:

- `GSE97362` and `GSE218186` use `process_meth = "funnorm"`.
- `GSE116300` uses `process_meth = "quantile"` and is also compared with the processed beta-values distributed through GEO.

This difference is important rather than incidental. Preprocessing can shift methylation distributions, so a model should not be assumed to transfer unchanged across preprocessing methods or platforms. In this tutorial, these differences are kept visible because the Kabuki analysis explicitly tests robustness across independent datasets and technologies.

## Probe filters

After normalization, the scripts apply three probe-level filters:

1. remove probes overlapping common SNPs with `dropLociWithSnps(maf = 0.01)`;
2. remove known cross-reactive probes with `maxprobes::dropXreactiveLoci()`;
3. keep only probes with detection P-value below 0.01 in every sample for the strict version of the workflow.

These filters reduce measurements that are more likely to be unreliable or difficult to interpret.


``` r
GRset <- minfi::dropLociWithSnps(GRset, maf = 0.01)
GRset <- maxprobes::dropXreactiveLoci(GRset)

detP <- detP[rownames(GRset), GRset$idat_name, drop = FALSE]
keep <- rowSums(detP < 0.01, na.rm = TRUE) == ncol(GRset)
print(table("CpGs kept" = keep))

GRset <- GRset[keep, ]
GRset$Dataset <- geo_ID
```

Requiring detection in every sample is conservative. It minimizes missing values but can remove many CpGs as cohort size grows. The tutorial helper can instead use a 95% detected-sample rule when a small amount of missingness is acceptable. Whichever rule is used should be fixed before downstream analysis and reported clearly.

## Wrapper function

For a new dataset, the wrapper collects these preprocessing choices in one function and returns a normalized `GenomicRatioSet` for downstream analysis. This makes it easier to see which settings were used for each cohort.


``` r
GRset <- read_idat_minfi(
  path_idat = path_idat,
  pheno = pheno,
  path_save = path_save,
  estimate_cell_counts = TRUE,
  dataset_name = geo_ID,
  process_meth = "funnorm"
)

plot_custom_theme <- function() {
  theme_classic() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 17, face = "bold", hjust = 0.5)
    )
}

saveRDS(GRset, file.path(path_save, paste0("GRset_", geo_ID, ".rds")))
```

## Cell composition

Whole-blood DNA contains a mixture of leukocyte populations, and each cell type has a characteristic methylation profile. Differences in cell composition between cases and controls can therefore create apparent methylation differences that are not directly related to the disorder.

The tutorial uses `FlowSorted.Blood.EPIC` to estimate major immune-cell fractions from the raw `RGChannelSet` [@salas2018blood]. These estimates are later included as covariates in the Williams differential methylation model. Cell deconvolution is one way to reduce confounding; it does not make the samples biologically identical.



``` r
RGset <- minfi::read.metharray.exp(base = path_idat) # Reads the IDATs.

library(FlowSorted.Blood.EPIC)

cell_estimates <- estimateCellCounts2(
  RGset, compositeCellType = "Blood",
  processMethod = "preprocessNoob",
  probeSelect = "IDOL",
  cellTypes = c("CD8T", "CD4T", "NK",
                "Bcell","Mono", "Neu")
)

as.data.frame(cell_estimates$prop) %>%
  rownames_to_column(var = "geo_accession") %>%
  mutate(
    "geo_accession" = sub("_.*", "", geo_accession)
  ) %>%
  readr::write_csv(file = paste0(path_save, "Tables/CellCounts_", geo_ID, ".csv"))
```

## PCA for batch effects

Before testing for differential methylation, the workflow uses PCA to inspect the largest sources of variation. For this QC view, variability is calculated from the 10,000 most variable autosomal CpGs among controls, and the resulting PCA is then inspected against slide, array position, sex, and case/control status. The purpose is to identify possible technical structure or unusual samples before feature discovery.


``` r
GRset$Slide <- stringr::str_split(GRset$idat_name, "_") %>%
  purrr::map_chr(~ .x[2]) %>%
  factor()

GRset$Array <- stringr::str_split(GRset$idat_name, "_") %>%
  purrr::map_chr(~ .x[3]) %>%
  factor()

annots <- if (exists("array_annots", mode = "function")) {
  array_annots(GRset)
} else {
  minfi::getAnnotation(GRset)
}
rowData(GRset) <- annots[rownames(GRset), ]

GRset_noXY <- GRset[!SummarizedExperiment::rowData(GRset)$chr %in% c("chrX", "chrY"), ]

variance_m <- apply(
  minfi::getM(GRset_noXY[, GRset_noXY$Case_Cont == "Control"]),
  1,
  stats::var
)
top_var <- sort(variance_m, decreasing = TRUE)[1:10000]

pca_res <- stats::prcomp(
  t(minfi::getBeta(GRset_noXY[names(top_var), ])),
  scale. = TRUE
)

var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)

pca_df <- as.data.frame(pca_res$x) %>%
  dplyr::select(PC1, PC2) %>%
  tibble::rownames_to_column("geo_accession") %>%
  dplyr::left_join(as.data.frame(SummarizedExperiment::colData(GRset_noXY)), by = "geo_accession")
```

The same PCA coordinates can then be coloured by slide, array, and case/control status to see whether one of these variables dominates the first principal components:


``` r
pca_panels <- lapply(c("Slide", "Array", "Case_Cont"), function(col_var) {
  ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2)) +
    ggplot2::stat_ellipse(
      ggplot2::aes(color = .data[[col_var]]),
      level = 0.95,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(ggplot2::aes(color = .data[[col_var]], shape = Sex), size = 2.5) +
    ggplot2::labs(
      title = paste("PCA plot -", col_var),
      x = paste0("PC1 (", var_explained[1], "% var)"),
      y = paste0("PC2 (", var_explained[2], "% var)")
    ) +
    plot_custom_theme()
})

ggplot2::ggsave(
  file.path(path_figures, paste0("PCA_", geo_ID, ".png")),
  ggpubr::ggarrange(plotlist = pca_panels, ncol = 3),
  width = 13,
  height = 5,
  bg = "white"
)
```

## PCA outliers and held-out validation

For the Williams example, the script uses the PC1/PC2 coordinates to flag samples outside the group-specific PCA distribution and also reserves five non-outlier samples per case/control group as hold-outs. Together, the reviewed outliers and the random hold-outs form the within-dataset validation set used later in the tutorial.


``` r
outliers_df <- detect_pca_outliers(pca_df, return_all = TRUE)

outliers_df %>%
  dplyr::filter(outside) %>%
  dplyr::select(geo_accession, Case_Cont, D2, outside) %>%
  openxlsx::write.xlsx(
    file.path(path_save, "Tables", paste0("PCA_outliers_", geo_ID, ".xlsx")),
    rowNames = FALSE
  )

set.seed(42)
outliers_df %>%
  dplyr::filter(!outside) %>%
  dplyr::group_by(Case_Cont) %>%
  dplyr::slice_sample(n = min(5, dplyr::n())) %>%
  dplyr::select(geo_accession, Case_Cont, D2, outside) %>%
  openxlsx::write.xlsx(
    file.path(path_save, "Tables", paste0("KeepOut_Validation_", geo_ID, ".xlsx")),
    rowNames = FALSE
  )
```

PCA outlier flags are prompts for review, not automatic exclusion rules. Before removing a sample, check its phenotype, array position, signal quality, sex prediction, and any other available technical information. This is especially important in rare-disease cohorts, where a biologically unusual patient may also be scientifically informative.

