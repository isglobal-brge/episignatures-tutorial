# Wrapper functions {#wrapper-functions}

This chapter collects the reusable R functions used throughout the tutorial. Keeping these functions in one place makes the practical chapters easier to read and gives readers a single reference for the code that supports preprocessing, model training, prediction, and visualization.

The functions shown here correspond to the helper script used by the tutorial:

```text
R/functions_tutorial_epis.R
```

In the practical chapters, that file is loaded with:


``` r
source(file.path("R", "functions_tutorial_epis.R"))
```

The code in this chapter is shown with `eval=FALSE`. This means that **rendering the book displays the functions but does not execute them**. Analyses should still use the sourceable `.R` file above.

Several small safeguards are included to make the wrappers easier to reuse outside the original analysis environment. In particular, output paths are passed explicitly, required model CpGs are checked before prediction, and plotting functions avoid machine-specific temporary paths.

## Project and package helpers

These small functions are used internally by the rest of the wrapper file. `.require_package()` gives a clear error if a required package is unavailable, while `ensure_project_dirs()` creates the output folders expected by the tutorial.


``` r
.require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Please install it first.", call. = FALSE)
  }
  invisible(TRUE)
}

ensure_project_dirs <- function(path_save) {
  dirs <- c(
    path_save,
    file.path(path_save, "IDATs"),
    file.path(path_save, "Tables"),
    file.path(path_save, "Figures"),
    file.path(path_save, "tmp_data")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}
```

## Reading and preprocessing methylation arrays

`read_idat_minfi()` performs the main preprocessing steps used for the array datasets in the tutorial. It reads raw IDAT files, applies the sample-level detection P-value filter, joins phenotype information, predicts sex for quality control, normalizes the methylation data, removes probes overlapping common SNPs and cross-reactive probes, and applies the probe-level detection P-value filter.

When requested, it also estimates blood-cell proportions and writes them to the tutorial output directory. The function returns a processed `GenomicRatioSet` for downstream analyses.

The quality-control thresholds used here reproduce the tutorial workflow and should not be interpreted as universal thresholds for every methylation study.


``` r
read_idat_minfi <- function(
  path_idat,
  pheno,
  path_save = NULL,
  estimate_cell_counts = TRUE,
  dataset_name = "",
  process_meth = "funnorm"
) {
  .require_package("minfi")
  .require_package("dplyr")
  .require_package("tibble")
  .require_package("S4Vectors")
  .require_package("SummarizedExperiment")
  .require_package("maxprobes")

  if (!dir.exists(path_idat)) {
    stop("IDAT directory does not exist: ", path_idat, call. = FALSE)
  }
  if (!"geo_accession" %in% colnames(pheno)) {
    stop("`pheno` must contain a `geo_accession` column.", call. = FALSE)
  }

  process_meth <- match.arg(process_meth, c("funnorm", "noob", "quantile"))

  message("Reading IDATs...")
  RGset <- minfi::read.metharray.exp(base = path_idat)

  message("Sample detection P-values...")
  detP <- minfi::detectionP(RGset)
  sample_detP <- colMeans(detP, na.rm = TRUE)
  keep_samples <- sample_detP < 0.05
  message("Keeping ", sum(keep_samples), " / ", length(keep_samples), " samples.")

  if (!any(keep_samples)) {
    stop("No samples passed the mean detection P-value < 0.05 filter.", call. = FALSE)
  }
  RGset <- RGset[, keep_samples]

  sample_meta <- data.frame(
    geo_accession = sub("_.*", "", colnames(RGset)),
    idat_name = colnames(RGset),
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(pheno, by = "geo_accession")

  rownames(sample_meta) <- sample_meta$geo_accession
  SummarizedExperiment::colData(RGset) <- S4Vectors::DataFrame(sample_meta)

  message("Predicting sex for QC...")
  MSet_qc <- minfi::preprocessNoob(RGset)
  MSet_qc <- minfi::mapToGenome(MSet_qc)

  predSex <- as.data.frame(minfi::getSex(MSet_qc)) |>
    tibble::rownames_to_column(var = "sample_id") |>
    dplyr::mutate(
      geo_accession = sub("_.*", "", sample_id),
      predSex = dplyr::case_when(
        predictedSex == "M" ~ "Male",
        predictedSex == "F" ~ "Female",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::select(geo_accession, predictedSex, predSex)

  sample_meta <- as.data.frame(SummarizedExperiment::colData(RGset)) |>
    tibble::rownames_to_column(var = ".row_id") |>
    dplyr::left_join(predSex, by = "geo_accession")
  rownames(sample_meta) <- sample_meta$.row_id
  sample_meta$.row_id <- NULL
  SummarizedExperiment::colData(RGset) <- S4Vectors::DataFrame(sample_meta)

  message("Processing the data...")
  if (process_meth == "funnorm") {
    message("Using functional normalization...")
    GRset <- minfi::preprocessFunnorm(RGset, sex = RGset$predictedSex)
  } else if (process_meth == "noob") {
    message("Using noob normalization...")
    GRset <- minfi::preprocessNoob(RGset)
  } else {
    message("Using quantile normalization...")
    GRset <- minfi::preprocessQuantile(RGset)
  }

  message("Dropping loci with SNPs...")
  GRset <- minfi::dropLociWithSnps(GRset, maf = 0.01)

  message("Dropping cross-reactive probes...")
  GRset <- maxprobes::dropXreactiveLoci(GRset)

  # For other array versions, consider platform-specific masks, e.g. sesame.

  message("Filtering probes with low detection P-values...")
  GRset <- GRset[intersect(rownames(GRset), rownames(detP)), ]
  detP_subset <- detP[rownames(GRset), GRset$idat_name, drop = FALSE]
  keep_probes <- rowSums(detP_subset < 0.01, na.rm = TRUE) == ncol(GRset)
  GRset <- GRset[keep_probes, ]

  GRset$Dataset <- dataset_name

  if (estimate_cell_counts) {
    .require_package("FlowSorted.Blood.EPIC")
    .require_package("readr")

    if (is.null(path_save)) {
      stop("`path_save` must be supplied when `estimate_cell_counts = TRUE`.", call. = FALSE)
    }
    ensure_project_dirs(path_save)

    message("Estimating immune cell composition...")
    cell_estimates <- FlowSorted.Blood.EPIC::estimateCellCounts2(
      RGset,
      compositeCellType = "Blood",
      processMethod = "preprocessNoob",
      probeSelect = "IDOL",
      cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
    )

    cell_counts <- as.data.frame(cell_estimates$prop) |>
      tibble::rownames_to_column(var = "geo_accession") |>
      dplyr::mutate(geo_accession = sub("_.*", "", geo_accession))

    readr::write_csv(
      cell_counts,
      file.path(path_save, "Tables", paste0("CellCounts_", dataset_name, ".csv"))
    )
  }

  message("Finished processing the data.")
  GRset
}
```

## Shared plotting theme

`plot_custom_theme()` provides a common `ggplot2` theme so figures have a consistent appearance across the tutorial.


``` r
plot_custom_theme <- function() {
  .require_package("ggplot2")
  ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right",
      axis.text = ggplot2::element_text(size = 13),
      axis.title = ggplot2::element_text(size = 15, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.title = ggplot2::element_text(size = 17, face = "bold", hjust = 0.5)
    )
}
```

## PCA-based outlier identification

`detect_pca_outliers()` calculates Mahalanobis distance using the first two principal components within each case/control group. In the Williams syndrome workflow, the resulting flag is used to identify atypical profiles that can be reserved for validation.

Importantly, an outlier flag is an exploratory statistical result. It should not by itself be interpreted as evidence that a biological sample is poor quality or should be discarded.


``` r
detect_pca_outliers <- function(
  pca_df,
  conf_level = 0.95,
  return_all = TRUE
) {
  .require_package("dplyr")

  required <- c("PC1", "PC2", "Case_Cont")
  missing <- setdiff(required, colnames(pca_df))
  if (length(missing) > 0) {
    stop("Missing PCA columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  chi_cutoff <- stats::qchisq(conf_level, df = 2)

  result <- pca_df |>
    dplyr::group_by(Case_Cont) |>
    dplyr::group_modify(function(.x, .y) {
      xy <- cbind(.x$PC1, .x$PC2)
      if (nrow(.x) < 3) {
        .x$D2 <- NA_real_
        .x$outside <- NA
        return(.x)
      }
      d2 <- tryCatch(
        stats::mahalanobis(
          xy,
          center = c(mean(.x$PC1), mean(.x$PC2)),
          cov = stats::cov(xy)
        ),
        error = function(e) rep(NA_real_, nrow(.x))
      )
      .x$D2 <- d2
      .x$outside <- d2 > chi_cutoff
      .x
    }) |>
    dplyr::ungroup()

  if (return_all) result else dplyr::filter(result, outside %in% TRUE)
}
```

## Preparing binary classification models

The following internal functions convert methylation values into the sample-by-feature format required by `caret`, create the binary case/control outcome, calculate class weights, perform repeated cross-validation for hyperparameter selection, and fit the final model.

The outcome is coded so that **Case** is the positive class. This is important when interpreting ROC-based tuning and the probability column returned during prediction.


``` r
.prepare_binary_df <- function(GRset_sbst, beta_M = "Beta") {
  .require_package("minfi")
  .require_package("SummarizedExperiment")

  beta_M <- match.arg(beta_M, c("Beta", "M"))

  if (!"Case_Cont" %in% colnames(SummarizedExperiment::colData(GRset_sbst))) {
    stop("`GRset_sbst` must contain a `Case_Cont` column.", call. = FALSE)
  }
  if (anyNA(GRset_sbst$Case_Cont)) {
    stop("`Case_Cont` contains missing values in the training set.", call. = FALSE)
  }
  if (!all(as.character(GRset_sbst$Case_Cont) %in% c("Case", "Control"))) {
    stop("`Case_Cont` must contain only 'Case' and 'Control'.", call. = FALSE)
  }

  message("Using assay: ", beta_M)
  X <- if (beta_M == "Beta") {
    t(minfi::getBeta(GRset_sbst))
  } else {
    t(minfi::getM(GRset_sbst))
  }

  df <- as.data.frame(X, check.names = FALSE)
  df$Diagnosis <- factor(
    as.character(GRset_sbst$Case_Cont),
    levels = c("Case", "Control")
  )
  df
}

.compute_class_weights <- function(y) {
  tbl <- table(y)
  if (!all(c("Case", "Control") %in% names(tbl))) {
    stop("Both Case and Control samples are required.", call. = FALSE)
  }
  n <- sum(tbl)
  w_case <- n / (2 * tbl[["Case"]])
  w_control <- n / (2 * tbl[["Control"]])
  w_case <- as.numeric(w_case / w_control)
  list(Case = w_case, Control = 1)
}

.train_binary_model <- function(
  GRset_sbst,
  method,
  tune_grid,
  repeats = 10,
  number = 5,
  sampling = NULL,
  pre_process = c("center", "scale"),
  n_cores = 5,
  beta_M = "Beta"
) {
  .require_package("caret")
  .require_package("doParallel")
  .require_package("foreach")

  message("*** Training ", method, " ***")
  train_df <- .prepare_binary_df(GRset_sbst, beta_M)

  class_weights <- .compute_class_weights(train_df$Diagnosis)
  sample_weights <- ifelse(
    train_df$Diagnosis == "Case",
    class_weights$Case,
    class_weights$Control
  )

  n_cores <- max(1L, as.integer(n_cores))
  if (n_cores > 1L) {
    cl <- parallel::makePSOCKcluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  tune_control <- caret::trainControl(
    method = "repeatedcv",
    number = number,
    repeats = repeats,
    classProbs = TRUE,
    sampling = sampling,
    summaryFunction = caret::twoClassSummary
  )

  extra_args <- list()
  if (method == "svmLinear2") {
    extra_args$class.weights <- class_weights
  } else if (method == "glmnet") {
    extra_args$weights <- sample_weights
  } else if (method == "ranger") {
    extra_args$case.weights <- sample_weights
  } else if (method == "nnet") {
    extra_args$trace <- FALSE
    extra_args$maxit <- 500
  }

  mod_cv <- do.call(caret::train, c(
    list(
      Diagnosis ~ .,
      data = train_df,
      method = method,
      trControl = tune_control,
      tuneGrid = tune_grid,
      preProcess = pre_process,
      metric = "ROC"
    ),
    extra_args
  ))

  best_grid <- mod_cv$bestTune

  final_control <- caret::trainControl(method = "none", classProbs = TRUE)
  mod_final <- do.call(caret::train, c(
    list(
      Diagnosis ~ .,
      data = train_df,
      method = method,
      trControl = final_control,
      tuneGrid = best_grid,
      preProcess = pre_process
    ),
    extra_args
  ))

  mod_final$class_weights <- class_weights
  mod_final$sample_weights <- sample_weights
  mod_final$beta_M <- beta_M
  mod_final$bestTune_cv <- best_grid
  mod_final$tuning_results <- mod_cv$results
  mod_final
}
```

## Machine-learning model wrappers

These functions provide a consistent interface for the classifiers evaluated in the Williams syndrome practical. They include a linear support vector machine, elastic-net regression, random forest through `ranger`, k-nearest neighbours, partial least squares, naive Bayes, gradient boosting, XGBoost, and a neural network.

The tutorial compares several algorithms to demonstrate classifier implementation. Strong performance in the Williams example reflects the strength of that particular methylation signal and should not be taken to mean that all algorithms will perform equally well for every episignature.


``` r
svm_training_f <- function(
  GRset_sbst,
  repeats = 10,
  number = 5,
  sampling = NULL,
  beta_M = "Beta",
  n_cores = 5
) {
  svm_grid <- expand.grid(
    cost = 10^seq(-4, 4, length.out = 40) # Cost for the linear SVM margin penalty.
  )
  .train_binary_model(
    GRset_sbst, "svmLinear2", svm_grid,
    repeats, number, sampling,
    c("center", "scale"), n_cores, beta_M
  )
}

elastic_net_training_f <- function(
  GRset_sbst, repeats = 10, number = 5, sampling = NULL,
  beta_M = "Beta", n_cores = 5
) {
  message("-- Elastic net training --")
  glmnet_grid <- expand.grid(
    alpha = seq(0, 1, by = 0.1),
    lambda = 10^seq(-5, 1, length.out = 50)
  )
  .train_binary_model(
    GRset_sbst, "glmnet", glmnet_grid,
    repeats, number, sampling, NULL, n_cores, beta_M
  )
}

ranger_training_f <- function(
  GRset_sbst, repeats = 10, number = 5, sampling = NULL,
  beta_M = "Beta", n_cores = 5
) {
  message("-- Ranger training --")
  p <- ncol(.prepare_binary_df(GRset_sbst, beta_M)) - 1
  ranger_grid <- expand.grid(
    mtry = unique(round(seq(sqrt(p), min(500, p), length.out = 8))),
    splitrule = "gini",
    min.node.size = c(1, 5, 10)
  )
  .train_binary_model(
    GRset_sbst, "ranger", ranger_grid,
    repeats, number, sampling, NULL, n_cores, beta_M
  )
}

knn_training_f <- function(
  GRset_sbst, repeats = 10, number = 5, sampling = NULL,
  beta_M = "Beta", n_cores = 5
) {
  message("-- KNN training --")
  train_df <- .prepare_binary_df(GRset_sbst, beta_M)
  n_samples <- nrow(train_df)
  n_minority <- min(table(train_df$Diagnosis))
  k_upper <- min(ceiling(1.5 * sqrt(n_samples)), floor(0.8 * n_minority))
  k_upper <- max(k_upper, 3)
  k_seq <- seq(3, k_upper, by = 2)
  if (length(k_seq) == 0) k_seq <- 3
  knn_grid <- expand.grid(kmax = k_seq, distance = 2, kernel = "optimal")
  .train_binary_model(
    GRset_sbst, "kknn", knn_grid,
    repeats, number, sampling,
    c("center", "scale"), n_cores, beta_M
  )
}

pls_training_f <- function(
  GRset_sbst, repeats = 10, number = 5, sampling = NULL,
  beta_M = "Beta", n_cores = 5
) {
  message("-- PLS training --")
  max_ncomp <- min(20, ncol(GRset_sbst) - 1, nrow(GRset_sbst))
  pls_grid <- expand.grid(ncomp = seq_len(max(1, max_ncomp)))
  .train_binary_model(
    GRset_sbst, "pls", pls_grid,
    repeats, number, sampling,
    c("center", "scale"), n_cores, beta_M
  )
}

nb_training_f <- function(
  GRset_sbst, repeats = 10, number = 5, sampling = NULL,
  beta_M = "Beta", n_cores = 5
) {
  message("-- Naive Bayes training --")
  nb_grid <- expand.grid(
    usekernel = c(TRUE, FALSE),
    fL = c(0, 0.5, 1),
    adjust = c(0.5, 1, 2, 3, 5)
  )
  .train_binary_model(
    GRset_sbst, "nb", nb_grid,
    repeats, number, sampling,
    c("center", "scale"), n_cores, beta_M
  )
}

gbm_training_f <- function(
  GRset_sbst, repeats = 10, number = 5, sampling = NULL,
  beta_M = "Beta", n_cores = 5
) {
  message("-- GBM training --")
  gbm_grid <- expand.grid(
    n.trees = c(100, 250, 500),
    interaction.depth = c(1, 2, 3),
    shrinkage = c(0.01, 0.05, 0.1),
    n.minobsinnode = 10
  )
  .train_binary_model(
    GRset_sbst, "gbm", gbm_grid,
    repeats, number, sampling, NULL, n_cores, beta_M
  )
}

xgb_training_f <- function(
  GRset_sbst, repeats = 10, number = 5, sampling = NULL,
  beta_M = "Beta", n_cores = 5
) {
  message("-- XGBoost training --")
  xgb_grid <- expand.grid(
    nrounds = c(100, 250),
    max_depth = c(2, 4, 6),
    eta = c(0.01, 0.05, 0.1),
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample = 0.8
  )
  .train_binary_model(
    GRset_sbst, "xgbTree", xgb_grid,
    repeats, number, sampling, NULL, n_cores, beta_M
  )
}

nnet_training_f <- function(
  GRset_sbst, repeats = 10, number = 5, sampling = NULL,
  beta_M = "Beta", n_cores = 5
) {
  message("-- Neural Network training --")
  nnet_grid <- expand.grid(
    size = c(1, 3, 5, 7, 10),
    decay = 10^seq(-5, -1, length.out = 5)
  )
  .train_binary_model(
    GRset_sbst, "nnet", nnet_grid,
    repeats, number, sampling,
    c("center", "scale"), n_cores, beta_M
  )
}
```

## Predicting new samples

External validation requires the new dataset to contain the same CpG predictors used to train the model. `model_predictor_names()` retrieves those CpGs, and `predict_model_on_grset()` checks that they are present before calculating class probabilities.

This check prevents a common error in which the feature set is silently reduced to the intersection of two datasets. A fitted classifier should not be applied after dropping predictors unless the model has been explicitly retrained or a validated imputation strategy has been defined.


``` r
model_predictor_names <- function(model) {
  predictors <- model$coefnames
  if (is.null(predictors) || length(predictors) == 0) {
    if (!is.null(model$trainingData)) {
      predictors <- setdiff(colnames(model$trainingData), ".outcome")
    }
  }
  if (is.null(predictors) || length(predictors) == 0) {
    stop("Could not determine model predictor names.", call. = FALSE)
  }
  predictors
}

predict_model_on_grset <- function(model, GRset, beta_M = NULL) {
  .require_package("minfi")
  .require_package("tibble")

  model_cpgs <- model_predictor_names(model)
  missing_cpgs <- setdiff(model_cpgs, rownames(GRset))

  if (length(missing_cpgs) > 0) {
    stop(
      "External data are missing ", length(missing_cpgs),
      " model CpG(s). Missing examples: ",
      paste(utils::head(missing_cpgs, 10), collapse = ", "),
      ". Retrain on common CpGs or define a validated imputation strategy.",
      call. = FALSE
    )
  }

  if (is.null(beta_M)) beta_M <- model$beta_M
  if (is.null(beta_M)) beta_M <- "Beta"
  beta_M <- match.arg(beta_M, c("Beta", "M"))

  X <- if (beta_M == "Beta") {
    t(minfi::getBeta(GRset[model_cpgs, , drop = FALSE]))
  } else {
    t(minfi::getM(GRset[model_cpgs, , drop = FALSE]))
  }

  newdata <- as.data.frame(X, check.names = FALSE)
  newdata <- newdata[, model_cpgs, drop = FALSE]

  probs <- stats::predict(model, newdata = newdata, type = "prob")
  pred <- stats::predict(model, newdata = newdata, type = "raw")

  out <- tibble::tibble(
    geo_accession = rownames(newdata),
    Prediction = as.character(pred)
  )
  if ("Case" %in% colnames(probs)) out$Case <- probs$Case
  if ("Control" %in% colnames(probs)) out$Control <- probs$Control
  out
}
```

## Plotting and heatmap functions

The final group of functions creates the classifier probability plots, hierarchical-clustering heatmaps, and categorical colour palettes used in the worked examples.

These functions are used for visualization. They do not determine the episignature or change classifier predictions.


``` r
svm_plt <- function(predictions_df, col_pal, shape_pal) {
  .require_package("ggplot2")
  .require_package("ggbeeswarm")

  required <- c("Status", "Case", "dataset_id")
  missing <- setdiff(required, colnames(predictions_df))
  if (length(missing) > 0) {
    stop("Missing plot columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  plot_df <- predictions_df
  status_levels <- names(col_pal)[names(col_pal) %in% unique(as.character(plot_df$Status))]
  plot_df$Status <- factor(plot_df$Status, levels = status_levels)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = Status, y = Case, fill = Status)) +
    ggplot2::geom_vline(
      xintercept = seq_along(status_levels),
      color = "gray", linewidth = 0.5, alpha = 0.75
    ) +
    ggplot2::geom_hline(
      yintercept = c(0, 0.25, 0.5, 0.75, 1),
      linetype = "dashed", color = "gray", linewidth = 0.3, alpha = 0.75
    ) +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(shape = dataset_id),
      size = 3, color = "black", varwidth = TRUE, width = 0.25, alpha = 0.75
    ) +
    ggplot2::labs(x = "Group", y = "Case probability") +
    ggplot2::scale_fill_manual(values = col_pal, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shape_pal) +
    plot_custom_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

hclust_plt <- function(
  betas,
  pheno,
  scale_betas = FALSE,
  target_names = FALSE,
  condition_cols,
  dataset_cols,
  sample_label_cols = NULL
) {
  .require_package("ComplexHeatmap")
  .require_package("circlize")
  .require_package("grid")

  mat <- as.data.frame(betas)
  meta <- as.data.frame(pheno)
  rownames(meta) <- meta$geo_accession
  meta <- meta[rownames(mat), , drop = FALSE]
  stopifnot(identical(rownames(meta), rownames(mat)))

  ha <- ComplexHeatmap::HeatmapAnnotation(
    Condition = meta$Status,
    Dataset = meta$dataset_id,
    Score = meta$Case,
    col = list(
      Condition = condition_cols,
      Dataset = dataset_cols,
      Score = circlize::colorRamp2(c(0, 0.5, 1), c("#D73027", "#F7F7F7", "#1A9850"))
    ),
    na_col = "grey90"
  )

  if (scale_betas) {
    message("Scaling betas")
    beta_col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2EC4B6", "#1D3557", "#FF4D8D"))
    mat <- t(scale(t(mat)))
  } else {
    message("Not scaling betas")
    beta_col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#8C510A", "#121111", "#01665E"))
    mat <- as.matrix(mat)
  }

  if (!isFALSE(target_names)) {
    all_names <- rownames(mat)
    all_names[!all_names %in% target_names] <- ""
    if (is.null(sample_label_cols)) {
      samples_col <- rep("black", nrow(mat))
    } else {
      samples_col <- unname(sample_label_cols[as.character(meta$Status)])
      samples_col[is.na(samples_col)] <- "black"
    }
  } else {
    all_names <- rep("", nrow(mat))
    samples_col <- rep("black", nrow(mat))
  }

  ht <- ComplexHeatmap::Heatmap(
    t(mat),
    name = "Betas",
    col = beta_col_fun,
    top_annotation = ha,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_labels = all_names,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 5, col = samples_col),
    heatmap_legend_param = list(
      title = "Betas",
      direction = "horizontal",
      legend_width = grid::unit(5, "cm")
    )
  )

  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "left"
  )
}

hclust_plt_williams <- function(
  betas,
  pheno,
  scale_betas = FALSE,
  target_names = FALSE,
  condition_cols
) {
  .require_package("ComplexHeatmap")
  .require_package("circlize")
  .require_package("grid")

  mat <- as.data.frame(t(betas))
  meta <- as.data.frame(pheno)
  rownames(meta) <- meta$geo_accession
  meta <- meta[rownames(mat), , drop = FALSE]
  stopifnot(identical(rownames(meta), rownames(mat)))

  ha <- ComplexHeatmap::HeatmapAnnotation(
    Condition = meta$Status,
    col = list(Condition = condition_cols)
  )

  if (scale_betas) {
    message("Scaling betas")
    beta_col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2EC4B6", "#1D3557", "#FF4D8D"))
    mat <- t(scale(t(mat)))
  } else {
    message("Not scaling betas")
    beta_col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#3939e0", "#121111", "#c01313"))
    mat <- as.matrix(mat)
  }

  ht <- ComplexHeatmap::Heatmap(
    t(mat),
    name = "Betas",
    col = beta_col_fun,
    top_annotation = ha,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    heatmap_legend_param = list(
      title = "Beta-value",
      direction = "horizontal",
      legend_width = grid::unit(5, "cm")
    )
  )

  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "left"
  )
}

make_brewer_palette <- function(f, palette = "Set2") {
  .require_package("RColorBrewer")
  f <- factor(f)
  n <- nlevels(f)
  if (n == 0) return(stats::setNames(character(0), character(0)))

  info <- RColorBrewer::brewer.pal.info
  if (!palette %in% rownames(info)) {
    stop("Unknown RColorBrewer palette: ", palette, call. = FALSE)
  }

  max_n <- info[palette, "maxcolors"]
  if (n <= max_n) {
    cols <- RColorBrewer::brewer.pal(max(3, n), palette)[seq_len(n)]
  } else {
    cols <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(max_n, palette)
    )(n)
  }
  stats::setNames(cols, levels(f))
}
```

## Using this file in the tutorial

For the Bookdown project, keep both files:

```text
R/functions_tutorial_epis.R
07-wrapper-functions.Rmd
```

The `.R` file is sourced by the analysis chapters. The `.Rmd` file is included in the book so readers can inspect the same supporting functions without leaving the tutorial.

If this chapter is listed in `_bookdown.yml`, for example:

```text
rmd_files:
  - index.Rmd
  - 01-primer.Rmd
  - 02-repository-map.Rmd
  - 03-setup-data.Rmd
  - 04-preprocessing-qc.Rmd
  - 05-williams-workflow.Rmd
  - 06-kabuki-workflow.Rmd
  - 07-wrapper-functions.Rmd
  - 10-references.Rmd
```

Bookdown will render it as the chapter following the Kabuki workflow. The filename controls its position in the book; the chapter title shown to readers comes from the first-level heading at the top of this file.
