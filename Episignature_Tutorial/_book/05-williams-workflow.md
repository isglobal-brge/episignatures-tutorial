# Williams Syndrome Workflow {#williams}

The Williams analysis is the de novo episignature-discovery example used in the manuscript. It starts with public blood methylation-array data, identifies case-control differentially methylated positions (DMPs), reduces those DMPs to a smaller CpG panel, trains classifiers, and then asks whether the selected pattern remains detectable in samples that were not used for discovery.

`GSE119778` is the discovery dataset from Kimura et al. [@kimura2020williams; @geoGSE119778]. `GSE66552`, from Strong et al., provides an independent external cohort containing typical controls, Williams syndrome samples, and 7q11.23 duplication samples [@strong2015williams; @geoGSE66552].

## Workflow overview

The main script is `scripts/Williams_analysis.R`. Age is not available in the `GSE119778` metadata used by this workflow, so age cannot be added to the differential methylation model. The model instead adjusts for the available sex information and estimated blood-cell proportions. This illustrates an important rule for small rare-disease datasets: use covariates that are available, scientifically justified, and supported by the sample size; do not add variables simply because they are commonly used elsewhere.

In the manuscript analysis, `GSE119778` contains 34 Williams syndrome samples and 34 controls before the discovery/validation split. Seven cases and seven controls are ultimately reserved for within-dataset validation, leaving 27 cases and 27 controls for discovery. The external `GSE66552` cohort contains 20 Williams syndrome samples, 15 typical controls, and 10 7q11.23 duplication samples.

The conceptual flow is:

1. Prepare `GSE119778` phenotype and IDAT data.
2. Normalize, filter, and save `GRset_GSE119778.rds`.
3. Merge estimated blood cell proportions.
4. Run PCA to inspect batch effects and select outliers/hold-outs.
5. Define a discovery set by excluding PCA outliers and validation hold-outs.
6. Run limma on M-values to find DMPs.
7. Rank and prune DMPs to create an episignature CpG panel.
8. Train classifiers on the discovery set.
9. Score discovery, internal validation, and `GSE66552` external samples.

## Prepare the discovery object

After the preprocessing steps in Chapter \@ref(preprocessing-qc), reload the normalized methylation object and attach the estimated blood-cell proportions that will be used as covariates in the DMP model.


``` r
geo_ID <- "GSE119778"
path_save <- file.path(project_root, "Williams", geo_ID)
path_figures <- file.path(path_save, "Figures")

GRset <- readRDS(file.path(path_save, paste0("GRset_", geo_ID, ".rds")))
GRset$Case_Cont <- factor(GRset$Case_Cont, levels = c("Control", "Case"))

cell_counts <- readr::read_csv(
  file.path(path_save, "Tables", paste0("CellCounts_", geo_ID, ".csv")),
  show_col_types = FALSE
)

colData(GRset) <- as.data.frame(SummarizedExperiment::colData(GRset)) %>%
  dplyr::left_join(cell_counts, by = "geo_accession") %>%
  tibble::column_to_rownames("geo_accession") %>%
  S4Vectors::DataFrame()
```

The later limma model includes the `CD8T`, `CD4T`, `NK`, `Bcell`, `Mono`, and `Neu` estimates. These columns therefore need to be attached to the sample metadata before the design matrix is created.

## Define discovery and validation samples

The PCA/QC step writes two tables:

- `PCA_outliers_GSE119778.xlsx`, containing samples flagged for review from the PCA;
- `KeepOut_Validation_GSE119778.xlsx`, containing the randomly selected non-outlier hold-outs.

Neither group is used for DMP discovery or CpG selection. This is important because the validation samples should remain independent of feature discovery [@ambroise2002selection].


``` r
outlier_smpl <- readxl::read_xlsx(
  file.path(path_save, "Tables", paste0("PCA_outliers_", geo_ID, ".xlsx"))
) %>%
  dplyr::pull(geo_accession)

validation_smpl <- readxl::read_xlsx(
  file.path(path_save, "Tables", paste0("KeepOut_Validation_", geo_ID, ".xlsx"))
) %>%
  dplyr::pull(geo_accession)

GRset_Discovery <- GRset[, !colnames(GRset) %in% c(outlier_smpl, validation_smpl)]
```

Remove sex chromosome probes before DMP testing:


``` r
annots <- minfi::getAnnotation(GRset_Discovery)
stopifnot(identical(rownames(annots), rownames(GRset_Discovery)))

GRset_Discovery <- GRset_Discovery[!annots$chr %in% c("chrX", "chrY"), ]
```

Removing X- and Y-chromosome probes reduces the chance that sex-associated methylation differences dominate the case-control DMP results. Reported/predicted sex is still included as a sample-level covariate in the model.

## Differential methylation with limma

The Williams script fits a limma model to M-values, with case/control status as the main term and sex plus estimated blood-cell fractions as covariates [@ritchie2015limma].


``` r
design <- stats::model.matrix(
  ~ Case_Cont + Sex + CD8T + CD4T + NK + Bcell + Mono + Neu,
  data = as.data.frame(SummarizedExperiment::colData(GRset_Discovery))
)

fit <- limma::lmFit(object = minfi::getM(GRset_Discovery), design = design)
fit <- limma::eBayes(fit = fit)
tt <- limma::topTable(fit, coef = "Case_ContCase", number = Inf)
```

The coefficient `Case_ContCase` means the contrast is **Case minus Control**, because `Case_Cont` was defined with `Control` as the reference level.

## Add interpretable effect sizes

The limma P-values are complemented with an effect size calculated from beta-values. For each CpG, the script calculates:

```text
deltaBeta = mean(beta in cases) - mean(beta in controls)
```

A positive `deltaBeta` means that the average beta-value is higher in Williams syndrome samples than in controls; a negative value means that it is lower. The M-values are used for the statistical model, whereas beta-values make the size and direction of the methylation difference easier to interpret [@du2010comparison].


``` r
deltaBeta <- as.data.frame(SummarizedExperiment::assay(GRset_Discovery, "Beta")) %>%
  tibble::rownames_to_column(var = "id") %>%
  tidyr::pivot_longer(cols = -id, names_to = "geo_accession", values_to = "Beta") %>%
  dplyr::left_join(
    as.data.frame(SummarizedExperiment::colData(GRset_Discovery)) %>%
      dplyr::select(geo_accession, Case_Cont),
    by = "geo_accession"
  ) %>%
  dplyr::group_by(id) %>%
  dplyr::reframe(
    deltaBeta = mean(Beta[Case_Cont == "Case"]) - mean(Beta[Case_Cont == "Control"])
  )

tt_annot <- tt %>%
  tibble::rownames_to_column(var = "id") %>%
  dplyr::left_join(deltaBeta, by = "id") %>%
  dplyr::select(id, P.Value, adj.P.Val, deltaBeta) %>%
  dplyr::left_join(
    as.data.frame(annots) %>%
      dplyr::select(Name, chr, pos, Relation_to_Island, UCSC_RefGene_Name),
    by = c("id" = "Name")
  )

openxlsx::write.xlsx(
  tt_annot,
  file.path(path_save, "Tables", paste0("DMPs_limma_", geo_ID, ".xlsx")),
  rowNames = FALSE
)
```

## Volcano plot

For the Williams example, the manuscript and tutorial use two filters:

- adjusted P-value (`adj.P.Val`) < 0.01;
- absolute methylation difference (`abs(deltaBeta)`) > 0.10.

These are analysis choices for this example, not universal episignature thresholds. In the manuscript analysis, they identify 1,681 hypermethylated and 1,082 hypomethylated DMPs before the additional feature-selection steps.


``` r
thresholds_DMPs <- list(adj.P.Val = 0.01, deltaBeta = 0.1)

plt_volcano <- tt_annot %>%
  dplyr::mutate(
    logP = -log10(adj.P.Val),
    sign = factor(
      dplyr::case_when(
        adj.P.Val < thresholds_DMPs[["adj.P.Val"]] &
          deltaBeta > thresholds_DMPs[["deltaBeta"]] ~ "HyperM",
        adj.P.Val < thresholds_DMPs[["adj.P.Val"]] &
          deltaBeta < -thresholds_DMPs[["deltaBeta"]] ~ "HypoM",
        TRUE ~ "NS"
      ),
      levels = c("HypoM", "NS", "HyperM")
    )
  ) %>%
  ggplot2::ggplot(ggplot2::aes(deltaBeta, logP, color = sign)) +
  ggplot2::geom_point(alpha = 0.2, size = 2) +
  ggplot2::geom_vline(xintercept = c(-0.1, 0.1), lty = 2, color = "gray") +
  ggplot2::geom_hline(yintercept = -log10(0.01), lty = 2, color = "gray") +
  ggplot2::scale_color_manual(values = c(NS = "gray", HyperM = "#c01313", HypoM = "#3939e0")) +
  ggplot2::labs(x = "Delta beta (Case - Control)", y = "-log10(FDR)", color = "") +
  plot_custom_theme()

ggplot2::ggsave(
  file.path(path_figures, paste0("Volcano_plt", geo_ID, ".png")),
  plt_volcano,
  width = 5,
  height = 4,
  bg = "white"
)
```

## Select episignature CpGs

The next steps reduce the significant DMPs to a smaller set of candidate classifier features. All of these decisions are made with the discovery samples only.

### Step 1: rank significant DMPs


``` r
tt_flt <- tt_annot %>%
  dplyr::filter(
    adj.P.Val < thresholds_DMPs[["adj.P.Val"]],
    abs(deltaBeta) > thresholds_DMPs[["deltaBeta"]]
  ) %>%
  dplyr::mutate(Pval_deltaBeta = -log10(adj.P.Val) * abs(deltaBeta)) %>%
  dplyr::arrange(dplyr::desc(Pval_deltaBeta)) %>%
  dplyr::slice_head(n = 1000)
```

The combined score prioritizes CpGs that have both strong statistical evidence and a large case-control methylation difference. It is a ranking device; it does not by itself show that an individual CpG has a known biological function.

### Step 2: remove redundant CpGs

Highly correlated CpGs carry overlapping information. Removing some of this redundancy reduces the number of predictors and helps prevent the final panel from being dominated by groups of CpGs with very similar behaviour, whether or not those CpGs are physically close in the genome.


``` r
beta_values <- minfi::getBeta(GRset_Discovery[tt_flt$id, ])

rm_CpG <- caret::findCorrelation(
  t(beta_values),
  cutoff = 0.9,
  verbose = TRUE,
  names = TRUE
)
```

### Step 3: keep high-importance CpGs

With `nonpara = TRUE`, `caret::filterVarImp()` evaluates each CpG individually using a two-class ROC-based measure of discrimination [@kuhn2008caret]. The script then keeps CpGs with importance above the median, corresponding to the more informative half of the remaining features.


``` r
var_Imp <- caret::filterVarImp(
  x = as.data.frame(t(beta_values[!rownames(beta_values) %in% rm_CpG, ])),
  y = GRset_Discovery$Case_Cont,
  nonpara = TRUE
) %>%
  tibble::rownames_to_column(var = "id") %>%
  dplyr::arrange(dplyr::desc(Case)) %>%
  dplyr::filter(Case > stats::quantile(Case, 0.5)) %>%
  dplyr::select(id, Case) %>%
  dplyr::rename(Importance = Case)

tt_flt_Imp <- tt_flt %>%
  dplyr::left_join(var_Imp, by = "id") %>%
  dplyr::filter(!is.na(Importance))

openxlsx::write.xlsx(
  tt_flt_Imp,
  file.path(path_save, "Tables", paste0("EpisignatureCpGs_", geo_ID, ".xlsx")),
  rowNames = FALSE
)
```

In the manuscript analysis, the filtering sequence reduces the top 1,000 ranked DMPs to 763 CpGs after correlation pruning and then to 374 CpGs after the variable-importance filter. These 374 CpGs form the candidate Williams episignature panel used for the later PCA and classifier training.

## Visualize the selected CpGs

The selected CpGs are first inspected with an unsupervised PCA. Because PCA does not use the case/control labels to create the axes, it provides a simple visual check of whether the selected methylation pattern separates the groups and where the held-out samples fall. This visualization supports the classifier analysis but is not a substitute for independent prediction.


``` r
pca_res <- stats::prcomp(
  t(minfi::getBeta(GRset[tt_flt_Imp$id, ])),
  scale. = TRUE,
  center = TRUE
)

var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)

pca_df <- as.data.frame(pca_res$x) %>%
  dplyr::select(PC1, PC2) %>%
  tibble::rownames_to_column(var = "geo_accession") %>%
  dplyr::left_join(as.data.frame(SummarizedExperiment::colData(GRset)), by = "geo_accession") %>%
  dplyr::mutate(
    Status = dplyr::case_when(
      !geo_accession %in% c(outlier_smpl, validation_smpl) & Case_Cont == "Case" ~ "Williams_D",
      !geo_accession %in% c(outlier_smpl, validation_smpl) & Case_Cont == "Control" ~ "Control_D",
      geo_accession %in% c(outlier_smpl, validation_smpl) & Case_Cont == "Case" ~ "Williams_V",
      geo_accession %in% c(outlier_smpl, validation_smpl) & Case_Cont == "Control" ~ "Control_V",
      TRUE ~ "Other"
    )
  )

ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = Status, shape = Sex)) +
  ggplot2::geom_point(size = 2.5, alpha = 0.8) +
  ggplot2::labs(
    title = "Williams episignature",
    x = paste0("PC1 (", var_explained[1], "% var)"),
    y = paste0("PC2 (", var_explained[2], "% var)")
  ) +
  plot_custom_theme()
```

## Train classifiers

The original Williams script compares several classifiers and trains the available `_training_f` wrappers using both beta-values and M-values. Depending on the installed packages, these include SVM, penalized/elastic-net models, random forest, KNN, PLS, naive Bayes, gradient boosting, XGBoost, and neural-network models.

For a first pass through the tutorial, the linear SVM is the simplest model to follow and is also consistent with the classifier family commonly used in episignature studies. The selected CpG panel is fixed before model training and before any validation samples are scored:


``` r
svm_mod <- svm_training_f(
  GRset_sbst = GRset_Discovery[tt_flt_Imp$id, ],
  repeats = 10,
  number = 5,
  sampling = NULL,
  beta_M = "Beta",
  n_cores = 5
)

saveRDS(svm_mod, file.path(path_save, paste0("svm_models_", geo_ID, ".rds")))
```

To reproduce the broader model sweep:


``` r
model_funs <- grep("_training_f$", ls(), value = TRUE)

mods_beta <- purrr::set_names(model_funs) %>%
  purrr::map(function(fun_name) {
    fun <- get(fun_name)
    res <- fun(
      GRset_sbst = GRset_Discovery[tt_flt_Imp$id, ],
      repeats = 10,
      number = 5,
      sampling = NULL,
      beta_M = "Beta"
    )
    res$algorithm_name <- fun_name
    res$beta_M <- "Beta"
    res
  })

mods_m <- purrr::set_names(model_funs) %>%
  purrr::map(function(fun_name) {
    fun <- get(fun_name)
    res <- fun(
      GRset_sbst = GRset_Discovery[tt_flt_Imp$id, ],
      repeats = 10,
      number = 5,
      sampling = NULL,
      beta_M = "M"
    )
    res$algorithm_name <- fun_name
    res$beta_M <- "M"
    res
  })

mods <- c(
  stats::setNames(mods_beta, paste0(names(mods_beta), "_Beta")),
  stats::setNames(mods_m, paste0(names(mods_m), "_M"))
)

saveRDS(mods, file.path(path_save, paste0("_classification_models_", geo_ID, ".rds")))
```

The larger model sweep is useful for comparing algorithms, but identical performance in this strong Williams example should not be interpreted as evidence that all algorithms are interchangeable. The manuscript analysis reports correct classification of the available Williams syndrome and control samples across the evaluated models in both internal and external validation; that result reflects the strength of this particular methylation signal and should not be generalized to other disorders.

## Score discovery and validation samples

For a list of models:


``` r
mods <- readRDS(file.path(path_save, paste0("_classification_models_", geo_ID, ".rds")))

GRset_Discovery <- GRset[, !colnames(GRset) %in% c(outlier_smpl, validation_smpl)]
GRset_Validation <- GRset[, !colnames(GRset) %in% colnames(GRset_Discovery)]

pred_discovery <- lapply(mods, function(spec_mod) {
  pred <- predict_model_on_grset(spec_mod, GRset_Discovery)
  pred %>%
    dplyr::left_join(as.data.frame(SummarizedExperiment::colData(GRset_Discovery)), by = "geo_accession") %>%
    dplyr::mutate(
      Status = dplyr::if_else(Case_Cont == "Case", "Williams_D", "Control_D"),
      Dataset = geo_ID,
      Algorithm = spec_mod$algorithm_name,
      beta_M = spec_mod$beta_M
    ) %>%
    dplyr::select(geo_accession, Status, Case, Dataset, Algorithm, beta_M)
})

pred_validation <- lapply(mods, function(spec_mod) {
  pred <- predict_model_on_grset(spec_mod, GRset_Validation)
  pred %>%
    dplyr::left_join(as.data.frame(SummarizedExperiment::colData(GRset_Validation)), by = "geo_accession") %>%
    dplyr::mutate(
      Status = dplyr::if_else(Case_Cont == "Case", "Williams_V", "Control_V"),
      Dataset = geo_ID,
      Algorithm = spec_mod$algorithm_name,
      beta_M = spec_mod$beta_M
    ) %>%
    dplyr::select(geo_accession, Status, Case, Dataset, Algorithm, beta_M)
})
```

## Prepare and score `GSE66552`

The independent `GSE66552` preparation script creates:

```text
Williams/GSE66552/GRset_epi_GSE66552.rds
```

This object is restricted to the Williams CpG panel selected from `GSE119778`. The external samples were not used to choose those CpGs or fit the classifiers.


``` r
geo_ID_external <- "GSE66552"
GRset_external <- readRDS(
  file.path(project_root, "Williams", geo_ID_external, paste0("GRset_epi_", geo_ID_external, ".rds"))
)

pred_external <- lapply(mods, function(spec_mod) {
  pred <- predict_model_on_grset(spec_mod, GRset_external)
  pred %>%
    dplyr::left_join(as.data.frame(SummarizedExperiment::colData(GRset_external)), by = "geo_accession") %>%
    dplyr::mutate(
      Status = dplyr::case_when(
        Diagnosis == "WS" ~ "Williams_E",
        Diagnosis == "TD control" ~ "Control_E",
        Diagnosis == "Dup7" ~ "Dup7_E",
        TRUE ~ NA_character_
      ),
      Dataset = geo_ID_external,
      Algorithm = spec_mod$algorithm_name,
      beta_M = spec_mod$beta_M
    ) %>%
    dplyr::select(geo_accession, Status, Case, Dataset, Algorithm, beta_M)
})
```

Combine all predictions:


``` r
predictions_df <- dplyr::bind_rows(
  dplyr::bind_rows(pred_discovery),
  dplyr::bind_rows(pred_validation),
  dplyr::bind_rows(pred_external)
) %>%
  dplyr::mutate(
    Status = factor(
      Status,
      levels = c("Control_D", "Williams_D", "Control_V", "Williams_V", "Control_E", "Williams_E", "Dup7_E")
    ),
    dataset_id = factor(Dataset, levels = c("GSE119778", "GSE66552"))
  )

openxlsx::write.xlsx(
  predictions_df,
  file.path(path_save, "Tables", paste0("Classification_predictions_", geo_ID, ".xlsx")),
  rowNames = FALSE
)
```

## Interpret the Williams outputs

When reviewing the final scores, look for the following pattern:

- discovery controls should score toward the control end and discovery Williams samples toward the case end;
- held-out samples should behave similarly to their expected group;
- independent Williams samples should retain the disease-associated signal;
- independent typical controls should remain control-like;
- 7q11.23 duplication samples provide an additional comparator and should not simply reproduce the Williams deletion profile.

Warning signs include a whole external dataset shifting together, much poorer performance in held-out samples than in discovery samples, scores that track a technical variable such as slide or array, or missing model CpGs in an external dataset.

The main lesson is that the Williams workflow is more than model fitting. It moves in order from discovery, to feature selection, to locked model training, and finally to independent testing. That separation is what allows the final result to address whether the methylation pattern generalizes beyond the samples that created it.



