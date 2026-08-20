# Review Notes and Things to Correct {#known-issues}

While reviewing the repository scripts, I found several places that are worth correcting or at least documenting before sharing the workflow broadly. None of these invalidate the whole project, but they are the kinds of details that can confuse new users or make reruns fragile.

## Hidden `path_save` dependency in the original helper

In `scripts/functions_TutorialEpis.R`, `read_idat_minfi()` writes cell counts using a global `path_save`:

```text
scripts/functions_TutorialEpis.R:99
```

That means the function can fail or write to the wrong folder if `path_save` is not defined in the calling environment. The bookdown helper fixes this by adding `path_save` as an explicit argument.

Recommended fix in the original helper:


```r
read_idat_minfi <- function(path_idat, pheno, path_save, estimate_cell_counts = TRUE, ...) {
  # ...
  readr::write_csv(
    cell_counts,
    file = file.path(path_save, "Tables", paste0("CellCounts_", dataset_name, ".csv"))
  )
}
```

## Kabuki deviation typo: `sd` should be `sd_beta`

In `scripts/Kabuki/Kabuki_analysis.R`, the first deviation block has:

```text
scripts/Kabuki/Kabuki_analysis.R:178
```

```r
z_score = beta_diff / sd
```

But `summary_betas` creates a column called `sd_beta`, not `sd`. Later in the same script, line 261 uses `sd_beta` correctly. The first block should be:


```r
diff_betas <- diff_betas %>%
  dplyr::mutate(
    beta_diff = Beta - mean_beta,
    z_score = beta_diff / sd_beta
  )
```

If the line with `sd` is executed, R will look for an object called `sd`; in many contexts that resolves to the base `sd()` function, which is not a numeric denominator.

## KNN training reads a missing `Diagnosis` column

In `scripts/functions_TutorialEpis.R`, `knn_training_f()` calculates minority class size using:

```text
scripts/functions_TutorialEpis.R:400
```

```r
n_minority <- min(table(GRset_sbst$Diagnosis))
```

Most `GenomicRatioSet` objects in this repository use `Case_Cont`, not `Diagnosis`, as the binary training label. The helper should calculate this from the prepared training data:


```r
train_df <- .prepare_binary_df(GRset_sbst, beta_M)
n_minority <- min(table(train_df$Diagnosis))
```

The new helper in this book uses that pattern.

## Comment mismatch in Williams CpG selection

In `scripts/Williams_analysis.R`, the variable-importance filter says:

```text
scripts/Williams_analysis.R:337
```

```r
Case > quantile(Case, 0.5) # Keep the top 25% most important CpGs.
```

`quantile(Case, 0.5)` keeps CpGs above the median, which is the top 50%, not the top 25%. Either update the comment or change the quantile to `0.75`.

## SVM grid comment mentions sigma, but the model is linear

In `scripts/functions_TutorialEpis.R`, the SVM grid comment says:

```text
scripts/functions_TutorialEpis.R:210
```

```r
cost = 10^seq(-4, 4, length.out = 40) # Sigma values ...
```

For `svmLinear2`, this grid tunes `cost`, not radial-kernel sigma. Recommended comment:

```r
# Cost values for the linear SVM margin penalty.
```

## Heatmap helper has hard-coded temp paths and a fragile regex

In `scripts/functions_TutorialEpis.R`, `hclust_plt()` writes temporary PDFs to:

```text
/PROJECTES/GENOMICS/aalegret/Epimutations_Episignatures/Datasets/tmp_plot.pdf
```

This makes the function machine-specific. It also contains:

```text
scripts/functions_TutorialEpis.R:644
```

```r
grepl("*.Array", names(samples_col))
```

The regex is fragile because `*` has special regex meaning. Use `grepl("Array", names(samples_col), fixed = TRUE)` or a stricter escaped pattern.

## Scripts are not clean first-run pipelines

Several dataset scripts begin by loading objects that are created later in the same file. For example:

```text
scripts/Williams_analysis.R:17
scripts/Williams_analysis.R:138
scripts/Kabuki/GSE97362/Kabuki_GSE97362_analysis.R:16
```

This is normal for exploratory R scripts run section-by-section, but confusing for readers who expect `source("script.R")` to work from a clean folder. For shareable workflows, split each dataset into:

1. `01_prepare_pheno.R`
2. `02_preprocess_idats.R`
3. `03_qc_pca.R`
4. `04_discover_signature.R`
5. `05_train_score.R`

Or use an explicit `targets`/`drake` pipeline.

## Hold-out sampling can fail in small groups

The Williams script samples five non-outlier samples per group:

```text
scripts/Williams_analysis.R:230
```

If a group has fewer than five available samples, this fails. Safer code:


```r
outliers_df %>%
  dplyr::filter(!outside) %>%
  dplyr::group_by(Case_Cont) %>%
  dplyr::slice_sample(n = min(5, dplyr::n()))
```

## Detection P-value thresholds should be documented

The helper and Williams script keep samples with mean detection P-value below 0.05:

```text
scripts/functions_TutorialEpis.R:24-25
scripts/Williams_analysis.R:60-61
```

This is a reasonable coarse filter, but community readers will ask whether sample failure was defined by mean detection P-value, fraction of failed probes, or another QC metric. The tutorial now describes the rule, but a methods section should state it explicitly.

## External prediction requires the full model CpG set

When scoring external data with `caret`, the new data must contain every predictor used during training. Intersecting CpGs across datasets is useful for visualization, but it is not enough for prediction if the intersection drops model CpGs.

Recommended guard:


```r
model_cpgs <- model_predictor_names(model)
missing <- setdiff(model_cpgs, rownames(GRset_external))
if (length(missing) > 0) {
  stop("External data are missing model CpGs. Retrain on common CpGs or define imputation.")
}
```

The helper `predict_model_on_grset()` already includes this guard.

## Reproducibility: set seeds around sampling and model training

The Williams script samples validation hold-outs and trains repeated cross-validation models. Set seeds before both:


```r
set.seed(42)
# hold-out sampling

set.seed(123)
# caret training
```

For exact reproducibility across parallel backends, consider using `trainControl(seeds = ...)`.

## Path portability

Most original scripts use absolute `/PROJECTES/...` paths. That is fine on the HPC filesystem, but not for a public tutorial. The new bookdown uses `params$project_root` and `file.path()` throughout so readers can adapt paths without search-and-replace.
