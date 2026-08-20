# Developing and Testing New Cases {#new-cases}

This chapter shows how the Williams and Kabuki examples can be adapted to a new syndrome, gene, validation cohort, or small case series. It is a research template rather than a clinical diagnostic pipeline, so each new project still needs its own design and validation plan.

## Before writing code

For a new episignature project, define the analysis design before selecting CpGs or training a classifier:

| Question | Why it matters |
|----------|----------------|
| What is the biological label? | The model can only learn the contrast represented by the labels. |
| What are the controls? | Controls should be appropriate for tissue, age, platform, and technical context whenever possible. |
| Are related disorders available? | Related disorders provide a stronger specificity test than healthy controls alone. |
| How many confirmed cases are available? | Rare-disease cohorts are often small, which increases overfitting risk. |
| Which tissue is being tested? | A blood episignature cannot automatically be assumed to exist in saliva, fibroblasts, or another tissue. |
| Which platform is being used? | 450K, EPIC, EPICv2, ONT, and PacBio differ in feature coverage and measurement properties. |

Before analysis, write down the discovery cohort, which samples will be held out, any independent cohorts, the planned DMP/CpG-selection rules, and how success will be judged. Keeping these decisions separate from the final test data reduces overfitting and selection bias [@ambroise2002selection].

## Minimal metadata contract

The template below uses a common phenotype schema. `Age` or other fields can be `NA` when they are genuinely unavailable, but the columns needed by the planned model must exist before that model is fitted:


``` r
required_columns <- c(
  "title",
  "geo_accession",
  "Sex",
  "Case_Cont",
  "Diagnosis",
  "Age",
  "Dataset"
)
```

This simple check catches missing columns and invalid binary labels before preprocessing begins. It does not verify that the biological labels themselves are correct, so the phenotype table should still be reviewed manually.


``` r
check_pheno <- function(pheno) {
  missing <- setdiff(required_columns, colnames(pheno))
  if (length(missing) > 0) {
    stop("Missing phenotype columns: ", paste(missing, collapse = ", "))
  }
  if (!all(pheno$Case_Cont %in% c("Case", "Control"))) {
    stop("Case_Cont must contain only Case and Control.")
  }
  invisible(TRUE)
}

check_pheno(pheno)
```

## Template: new discovery dataset

The following skeleton can be adapted to a syndrome-specific dataset. The main idea is to keep the file paths and metadata mapping explicit so that another reader can see exactly what was processed.


``` r
library(GEOquery)
library(minfi)
library(tidyverse)
library(limma)
library(caret)
library(doParallel)

source("scripts/bookdown_episignatures_tutorial/R/functions_tutorial_epis.R")

syndrome <- "MySyndrome"
geo_ID <- "GSEXXXXXX"
project_root <- "path/to/project/scripts"
path_save <- file.path(project_root, syndrome, geo_ID)
path_idat <- file.path(path_save, "IDATs")
path_figures <- file.path(path_save, "Figures")
ensure_project_dirs(path_save)

# 1. Metadata
gse <- GEOquery::getGEO(GEO = geo_ID, destdir = path_idat)[[1]]
raw_pheno <- Biobase::pData(gse)

pheno <- raw_pheno %>%
  dplyr::mutate(
    Sex = NA_character_,       # map from dataset-specific metadata
    Diagnosis = NA_character_, # map from dataset-specific metadata
    Case_Cont = NA_character_, # map from dataset-specific metadata
    Age = NA_real_,
    Dataset = geo_ID
  ) %>%
  dplyr::select(title, geo_accession, Sex, Case_Cont, Diagnosis, Age, Dataset)

check_pheno(pheno)
readr::write_csv(pheno, file.path(path_save, paste0("pheno_", geo_ID, ".csv")))

# 2. IDAT preprocessing
GRset <- read_idat_minfi(
  path_idat = path_idat,
  pheno = pheno,
  path_save = path_save,
  estimate_cell_counts = TRUE,
  dataset_name = geo_ID,
  process_meth = "funnorm"
)
saveRDS(GRset, file.path(path_save, paste0("GRset_", geo_ID, ".rds")))
```

## Template: DMP discovery


``` r
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

# After PCA review, load approved exclusions.
outlier_smpl <- readxl::read_xlsx(
  file.path(path_save, "Tables", paste0("PCA_outliers_", geo_ID, ".xlsx"))
) %>%
  dplyr::pull(geo_accession)

validation_smpl <- readxl::read_xlsx(
  file.path(path_save, "Tables", paste0("KeepOut_Validation_", geo_ID, ".xlsx"))
) %>%
  dplyr::pull(geo_accession)

GRset_Discovery <- GRset[, !colnames(GRset) %in% c(outlier_smpl, validation_smpl)]

annots <- minfi::getAnnotation(GRset_Discovery)
GRset_Discovery <- GRset_Discovery[!annots$chr %in% c("chrX", "chrY"), ]

design <- stats::model.matrix(
  ~ Case_Cont + Sex + CD8T + CD4T + NK + Bcell + Mono + Neu,
  data = as.data.frame(SummarizedExperiment::colData(GRset_Discovery))
)

fit <- limma::lmFit(minfi::getM(GRset_Discovery), design)
fit <- limma::eBayes(fit)
tt <- limma::topTable(fit, coef = "Case_ContCase", number = Inf)
```

If cell estimates are unavailable, do not leave those terms in the model. Use only covariates that are available and scientifically justified. Also avoid adding variables that contain no information, such as an age column that is entirely `NA`.


``` r
design <- stats::model.matrix(
  ~ Case_Cont + Sex + Age,
  data = as.data.frame(SummarizedExperiment::colData(GRset_Discovery))
)
```

## Template: CpG selection

The Williams procedure can be used as a worked starting point for DMP filtering and feature reduction. It should not be treated as a universal recipe for every disorder:


``` r
thresholds_DMPs <- list(adj.P.Val = 0.01, deltaBeta = 0.1)

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
  tibble::rownames_to_column("id") %>%
  dplyr::left_join(deltaBeta, by = "id")

tt_flt <- tt_annot %>%
  dplyr::filter(adj.P.Val < 0.01, abs(deltaBeta) > 0.1) %>%
  dplyr::mutate(score = -log10(adj.P.Val) * abs(deltaBeta)) %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::slice_head(n = 1000)

beta_values <- minfi::getBeta(GRset_Discovery[tt_flt$id, ])
rm_CpG <- caret::findCorrelation(t(beta_values), cutoff = 0.9, names = TRUE)

var_Imp <- caret::filterVarImp(
  x = as.data.frame(t(beta_values[!rownames(beta_values) %in% rm_CpG, ])),
  y = GRset_Discovery$Case_Cont,
  nonpara = TRUE
) %>%
  tibble::rownames_to_column("id") %>%
  dplyr::arrange(dplyr::desc(Case)) %>%
  dplyr::filter(Case > stats::quantile(Case, 0.5)) %>%
  dplyr::select(id, Importance = Case)

episignature_cpgs <- tt_flt %>%
  dplyr::inner_join(var_Imp, by = "id")
```

The thresholds in this example (`FDR < 0.01`, `|deltaBeta| > 0.10`, top 1,000 DMPs, correlation cutoff 0.90, and the upper half of the variable-importance distribution) reproduce the logic used in the Williams tutorial. They are not universal episignature parameters. A weaker methylation signal or a very small cohort may require a different strategy.

Whatever strategy is chosen, DMP testing, ranking, correlation filtering, and variable-importance selection must use the discovery/training data only. Validation samples should not influence which CpGs enter the model [@ambroise2002selection].

## Template: train and score


``` r
model <- svm_training_f(
  GRset_sbst = GRset_Discovery[episignature_cpgs$id, ],
  repeats = 10,
  number = 5,
  beta_M = "Beta",
  n_cores = 5
)

saveRDS(model, file.path(path_save, paste0("svm_model_", geo_ID, ".rds")))

scores_discovery <- predict_model_on_grset(model, GRset_Discovery)
```

The `predict_model_on_grset()` helper refuses to score a dataset when required model CpGs are missing. This is deliberate: silently dropping predictors changes the model being applied.

If a new platform does not measure the full CpG set, define the solution **before** looking at the external labels. For example, a new model can be trained on a pre-specified shared CpG set using the original training data only. Any imputation strategy should likewise be defined and validated in advance. The external test labels must not be used to choose the features or tune the replacement model.

## Template: external validation


``` r
external_geo <- "GSEYYYYYY"
path_external <- file.path(project_root, syndrome, external_geo)
GRset_external <- readRDS(file.path(path_external, paste0("GRset_", external_geo, ".rds")))

model <- readRDS(file.path(path_save, paste0("svm_model_", geo_ID, ".rds")))
scores_external <- predict_model_on_grset(model, GRset_external) %>%
  dplyr::left_join(as.data.frame(SummarizedExperiment::colData(GRset_external)), by = "geo_accession") %>%
  dplyr::mutate(dataset_id = external_geo)

readr::write_csv(
  scores_external,
  file.path(path_external, "Tables", paste0("SVM_scores_", external_geo, ".csv"))
)
```

Evaluate performance only after preprocessing, predictors, model parameters, and the scoring rule have been fixed. Useful outputs include:

- score distributions by diagnosis;
- ROC/AUC when a binary comparison is appropriate;
- sensitivity and specificity at a pre-specified threshold;
- a confusion matrix;
- plots separated by dataset and platform;
- PCA or heatmap views restricted to the episignature CpGs;
- specificity testing against related disorders when suitable samples are available.

For small rare-disease cohorts, report the individual sample behaviour as well as summary metrics; a single misclassified patient can be scientifically important even when an aggregate metric looks high.

## Adding a single new sample

For a new research sample, use the closest validated preprocessing and scoring workflow:

1. Place the IDAT pair in the appropriate dataset folder.
2. Create the required phenotype row.
3. Process the sample with suitable reference samples when the normalization method requires cohort context; single-sample processing can behave differently.
4. Review detection P-values, sex prediction when applicable, and availability of every required CpG.
5. Apply the locked model without changing the predictors or threshold in response to that sample.
6. Interpret the result against validated controls, relevant disease comparators, and known platform effects.

A negative episignature score does not by itself exclude pathogenicity. Variant-specific effects, mosaicism, biological heterogeneity, tissue dependence, or limited signature sensitivity can all contribute to a negative result [@kerkhof2024reporting].

## Reporting checklist

A community-facing episignature report should make the analysis reproducible and the limitations visible. Include:

- accession numbers or a data-availability statement;
- sample counts by group, sex, age, tissue, and platform when available;
- preprocessing and normalization method;
- sample- and probe-QC rules;
- outlier criteria and any manual exclusions;
- covariates used for differential methylation;
- the complete CpG-selection procedure;
- model type, tuning grid, cross-validation design, and random seed;
- the exact predictor CpG list and decision threshold;
- held-out and independent validation results;
- cross-disease specificity results when available;
- limitations, especially small cohort size, biological heterogeneity, tissue dependence, and platform transferability.

