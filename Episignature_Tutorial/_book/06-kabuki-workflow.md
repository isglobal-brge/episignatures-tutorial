# Kabuki Workflow {#kabuki}

The Kabuki workflow is the existing-signature validation and cross-platform example used in the manuscript. Unlike the Williams chapter, it does not discover a new CpG panel with limma. Instead, it starts from an established Kabuki methylation signature and tests whether the same disease-associated signal remains detectable across independent array cohorts and long-read methylation data [@arefeshghi2017kabuki; @arefeshghi2020episignatures].

The CpG list used by the manuscript contains 150 published sites. One site (`cg22948959`) is not available in all datasets, so the cross-dataset analysis uses the shared set available across the compared data. This is why CpG availability must be checked before a trained model is applied.

The core scripts are:

- `scripts/Kabuki/GSE97362/Kabuki_GSE97362_analysis.R`
- `scripts/Kabuki/GSE116300/Kabuki_GSE116300_analysis.R`
- `scripts/Kabuki/GSE218186/Kabuki_GSE218186_analysis.R`
- `scripts/Kabuki/PMID41225292/LR_analysis.R`
- `scripts/Kabuki/Kabuki_analysis.R`

A central lesson is that **successful classification is not the same as technical equivalence**. An array-trained classifier can preserve a strong Kabuki disease signal even when the numerical methylation values from arrays, ONT, and PacBio have different distributions. PCA, clustering, and the later deviation analysis are therefore used to examine platform effects separately from disease classification [@hildonen2026longread].

## Kabuki data sources

| Source | Role in workflow | Notes |
|--------|------------------|-------|
| `GSE97362` | array training cohort and within-cohort testing | Kabuki/CHARGE methylation study with confirmed KMT2D training labels [@butcher2017chargekabuki; @geoGSE97362] |
| `GSE116300` | independent array cohort | Kabuki cohort from Sobreira et al.; quantile normalization is used in this script [@sobreira2017kabuki; @geoGSE116300] |
| `GSE218186` | independent array cohort | reported sex is unavailable in the metadata used here, so predicted sex is stored [@geoGSE218186] |
| `PMID41225292` | ONT, PacBio, and matched array methylation data | methylation matrix from Hildonen et al., converted to a `GenomicRatioSet` [@hildonen2026longread] |

The fixed Kabuki CpG list used by the tutorial is stored at:

```text
scripts/Kabuki/PMID41225292/episignature_PMID41225292.xlsx
```

## Prepare the `GSE97362` training cohort

The `GSE97362` script first harmonizes the phenotype metadata from the original Kabuki/CHARGE study [@butcher2017chargekabuki].


``` r
geo_ID <- "GSE97362"
path_save <- file.path(project_root, "Kabuki", geo_ID)
path_idat <- file.path(path_save, "IDATs")
path_figures <- file.path(path_save, "Figures")
ensure_project_dirs(path_save)

gse <- GEOquery::getGEO(GEO = geo_ID, destdir = path_idat)[[1]]

pheno <- Biobase::pData(gse) %>%
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

readr::write_csv(pheno, file.path(path_save, paste0("pheno_", geo_ID, ".csv")))
```

Then it preprocesses IDATs:


``` r
GRset <- read_idat_minfi(
  path_idat = path_idat,
  pheno = readr::read_csv(file.path(path_save, paste0("pheno_", geo_ID, ".csv")), show_col_types = FALSE),
  path_save = path_save,
  estimate_cell_counts = FALSE,
  dataset_name = geo_ID,
  process_meth = "funnorm"
)

saveRDS(GRset, file.path(path_save, paste0("GRset_", geo_ID, ".rds")))
```

## Subset to published Kabuki CpGs


``` r
episignature_cpgs <- openxlsx::read.xlsx(
  file.path(project_root, "Kabuki", "PMID41225292", "episignature_PMID41225292.xlsx")
)

GRset_epis <- GRset[rownames(GRset) %in% episignature_cpgs$id, ]
```

The harmonized `colData` keeps the detailed publication labels in `Diagnosis` and `SampleType`, while also creating the binary `Case_Cont` field needed by the classifier. The detailed labels remain important because not every sample in the dataset is part of the training set.


``` r
colData(GRset_epis) <- as.data.frame(SummarizedExperiment::colData(GRset_epis)) %>%
  dplyr::mutate(
    Case_Cont = factor(
      dplyr::case_when(
        Diagnosis == "Control" ~ "Control",
        TRUE ~ "Case"
      ),
      levels = c("Case", "Control")
    )
  ) %>%
  dplyr::select(title, geo_accession, Case_Cont, Sex, Age, Diagnosis, Dataset, SampleType) %>%
  S4Vectors::DataFrame()

saveRDS(GRset_epis, file.path(path_save, paste0("GRset_epis_", geo_ID, ".rds")))
```

## Train the Kabuki SVM

The SVM is trained only on the molecularly confirmed KMT2D loss-of-function discovery subset and its designated controls. Other Kabuki-related, CHARGE/CHD7, and variant samples are kept outside this training subset so they can be evaluated afterward.


``` r
GRset_train <- GRset_epis[
  ,
  GRset_epis$SampleType %in% c(
    "Control for KMT2D LOF discovery cohort",
    "KMT2D LOF discovery cohort",
    "Control for validation cohort"
  )
]

GRset_train$Case_Cont <- factor(
  dplyr::case_when(
    GRset_train$SampleType == "Control for KMT2D LOF discovery cohort" ~ "Control",
    GRset_train$SampleType == "KMT2D LOF discovery cohort" ~ "Case",
    GRset_train$SampleType == "Control for validation cohort" ~ "Control",
    TRUE ~ NA_character_
  ),
  levels = c("Case", "Control")
)

svm_mod <- svm_training_f(GRset_sbst = GRset_train, beta_M = "Beta")

saveRDS(svm_mod, file.path(path_save, paste0("svm_model_", geo_ID, ".rds")))
```

The model returns a `Case` probability score. Higher values mean that the sample is more similar to the Kabuki class learned from the training data; lower values are more control-like. As in the Williams example, this is a model-derived score rather than a universal clinical cutoff.

## Score all `GSE97362` samples


``` r
svm_mod <- readRDS(file.path(path_save, paste0("svm_model_", geo_ID, ".rds")))

predictions_df <- predict_model_on_grset(svm_mod, GRset_epis) %>%
  dplyr::left_join(as.data.frame(SummarizedExperiment::colData(GRset_epis)), by = "geo_accession") %>%
  dplyr::mutate(
    DiscVal = factor(
      dplyr::case_when(
        SampleType %in% c("Control for KMT2D LOF discovery cohort", "Control for validation cohort") ~ "Train Control",
        SampleType == "KMT2D LOF discovery cohort" ~ "Train Case",
        TRUE ~ Diagnosis
      ),
      levels = c("Train Control", "Train Case", "Control", "KMT2D variant", "KDM6A variant", "CHARGE", "CHD7 variant")
    )
  )

openxlsx::write.xlsx(
  predictions_df,
  file.path(path_save, "Tables", paste0("SVM_predictions_", geo_ID, ".xlsx")),
  rowNames = FALSE
)
```

CHARGE/CHD7 samples provide a useful cross-disease specificity check because CHARGE and Kabuki can show overlapping developmental phenotypes but have distinct methylation signatures [@butcher2017chargekabuki]. A classifier that only separates Kabuki from healthy controls is a weaker test than one that is also challenged with a related disorder.

## Prepare `GSE116300`

`GSE116300`, an independent Kabuki cohort, is prepared with the same general object structure but uses quantile normalization in this script [@sobreira2017kabuki]:


``` r
geo_ID <- "GSE116300"
path_save <- file.path(project_root, "Kabuki", geo_ID)
path_idat <- file.path(path_save, "IDATs")
ensure_project_dirs(path_save)

GRset <- read_idat_minfi(
  path_idat = path_idat,
  pheno = readr::read_csv(file.path(path_save, paste0("pheno_", geo_ID, ".csv")), show_col_types = FALSE),
  path_save = path_save,
  estimate_cell_counts = FALSE,
  dataset_name = geo_ID,
  process_meth = "quantile"
)

GRset_epis <- GRset[rownames(GRset) %in% episignature_cpgs$id, ]

colData(GRset_epis) <- as.data.frame(SummarizedExperiment::colData(GRset_epis)) %>%
  dplyr::mutate(
    Diagnosis = mutation,
    Age = NA_real_,
    Sex = predSex
  ) %>%
  dplyr::select(title, geo_accession, Case_Cont, Sex, Age, Diagnosis, Dataset) %>%
  S4Vectors::DataFrame()

saveRDS(GRset_epis, file.path(path_save, paste0("GRset_epis_", geo_ID, ".rds")))
```

The script also compares beta-values generated from the IDAT workflow with the processed beta-values distributed through GEO. PCA provides a simple visual check of whether the two data sources show similar broad structure. This comparison is descriptive; it does not imply that the two preprocessing procedures produce identical methylation values.

## Prepare `GSE218186`

`GSE218186` provides another independent array-based testing cohort. It is not used to define the Kabuki CpG panel or to fit the SVM.


``` r
geo_ID <- "GSE218186"
path_save <- file.path(project_root, "Kabuki", geo_ID)
path_idat <- file.path(path_save, "IDATs")
ensure_project_dirs(path_save)

GRset <- read_idat_minfi(
  path_idat = path_idat,
  pheno = readr::read_csv(file.path(path_save, paste0("pheno_", geo_ID, ".csv")), show_col_types = FALSE),
  path_save = path_save,
  estimate_cell_counts = FALSE,
  dataset_name = geo_ID,
  process_meth = "funnorm"
)

GRset_epis <- GRset[rownames(GRset) %in% episignature_cpgs$id, ]

colData(GRset_epis) <- as.data.frame(SummarizedExperiment::colData(GRset_epis)) %>%
  dplyr::mutate(
    Diagnosis = dplyr::case_when(Case_Cont == "Case" ~ "Kabuki", TRUE ~ "Control"),
    Age = NA_real_,
    Sex = predSex
  ) %>%
  dplyr::select(title, geo_accession, Case_Cont, Sex, Age, Diagnosis, Dataset) %>%
  S4Vectors::DataFrame()

saveRDS(GRset_epis, file.path(path_save, paste0("GRset_epis_", geo_ID, ".rds")))
```

## Convert long-read beta workbook to a `GenomicRatioSet`

The long-read script reads the methylation workbook from Hildonen et al. [@hildonen2026longread]. In this file, a few rows contain sample metadata (`Platform`, `Affected`, and `Cohort`) and the remaining rows contain CpG methylation values. The code separates those two parts before creating the `GenomicRatioSet`.


``` r
path_save <- file.path(project_root, "Kabuki", "PMID41225292")

table_metadata <- readxl::read_excel(
  file.path(path_save, "Kabuki_PB_ONT_PMID41225292.xlsx")
)

pheno_lr <- table_metadata %>%
  dplyr::filter(ID %in% c("Platform", "Affected", "Cohort")) %>%
  tidyr::pivot_longer(cols = -ID, names_to = "geo_accession", values_to = "value") %>%
  tidyr::pivot_wider(names_from = ID, values_from = value) %>%
  dplyr::mutate(
    title = sub("_.*", "", geo_accession),
    Case_Cont = dplyr::case_when(Affected == "Affected" ~ "Case", TRUE ~ "Control"),
    Diagnosis = dplyr::case_when(
      Platform == "Array" & Case_Cont == "Case" ~ "Kabuki",
      Platform == "Array" & Case_Cont == "Control" ~ "Control",
      Case_Cont == "Case" ~ paste0("Kabuki_", Platform),
      Case_Cont == "Control" ~ paste0("Control_", Platform)
    ),
    Dataset = Platform,
    Sex = NA_character_,
    Age = NA_real_
  ) %>%
  dplyr::select(title, geo_accession, Case_Cont, Sex, Age, Diagnosis, Dataset) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("geo_accession") %>%
  S4Vectors::DataFrame()

betas_lr <- table_metadata %>%
  dplyr::filter(!ID %in% c("Platform", "Affected", "Cohort")) %>%
  tibble::column_to_rownames(var = "ID") %>%
  as.matrix()

GRset_lr <- minfi::makeGenomicRatioSetFromMatrix(
  betas_lr,
  pData = pheno_lr,
  what = "Beta"
)

saveRDS(GRset_lr, file.path(path_save, "GRset_epis_PMID41225292.rds"))
```

This object contains only the CpGs available in the workbook. A trained classifier can only be applied if all of its required predictor CpGs are present. The later guard performs this check explicitly rather than silently scoring an incomplete feature set.

## Combine Kabuki datasets

`Kabuki_analysis.R` loads the prepared `GRset_epis_*` objects and identifies the CpGs shared across them. This common set is useful for direct cross-dataset visualization, but prediction still requires the complete predictor set expected by the trained model.


``` r
path_save <- file.path(project_root, "Kabuki")
svm_mod <- readRDS(file.path(path_save, "GSE97362", "svm_model_GSE97362.rds"))

list_files <- list.files(
  path = path_save,
  pattern = "GRset_epis_",
  full.names = TRUE,
  recursive = TRUE
)

list_GRsets <- lapply(list_files, function(x) {
  GRset <- readRDS(x)
  if ("CN" %in% SummarizedExperiment::assayNames(GRset)) {
    SummarizedExperiment::assay(GRset, "CN") <- NULL
  }
  GRset
})

common_cpgs <- Reduce(intersect, lapply(list_GRsets, rownames))

model_cpgs <- model_predictor_names(svm_mod)
missing_model_cpgs <- setdiff(model_cpgs, common_cpgs)
if (length(missing_model_cpgs) > 0) {
  stop("The combined Kabuki object is missing model CpGs. Do not score until this is resolved.")
}

GRset_common <- Reduce(cbind, lapply(list_GRsets, function(x) x[common_cpgs, ]))
```

The explicit check prevents accidental scoring with missing features. If the combined data do not contain every CpG used to train the SVM, the code stops. In a real analysis, the solution must be defined before scoring—for example, training a new model on a pre-specified shared CpG set using the training data only. Missing predictors should not simply be ignored.

## Assign plot statuses

For plotting and interpretation, the script separates confirmed training samples from independent array samples, ONT/PacBio samples, and CHARGE/CHD7 comparators. These labels do not change the model; they make it possible to see whether prediction behaviour differs by disease group, dataset, or platform.


``` r
GRset_common$Status <- factor(
  dplyr::case_when(
    colnames(GRset_common) %in% rownames(svm_mod$trainingData) & GRset_common$Case_Cont == "Case" ~ "Train Kabuki",
    colnames(GRset_common) %in% rownames(svm_mod$trainingData) & GRset_common$Case_Cont == "Control" ~ "Train Control",
    GRset_common$Diagnosis %in% c("CHARGE", "CHD7 variant") ~ "CHARGE",
    GRset_common$Dataset == "ONT_A" & GRset_common$Case_Cont == "Case" ~ "Kabuki_ONT_A",
    GRset_common$Dataset == "ONT_WGS" & GRset_common$Case_Cont == "Case" ~ "Kabuki_ONT_WGS",
    GRset_common$Dataset == "PacBio" & GRset_common$Case_Cont == "Case" ~ "Kabuki_PacBio",
    GRset_common$Dataset == "ONT_A" & GRset_common$Case_Cont == "Control" ~ "Control_ONT_A",
    GRset_common$Dataset == "ONT_WGS" & GRset_common$Case_Cont == "Control" ~ "Control_ONT_WGS",
    GRset_common$Dataset == "PacBio" & GRset_common$Case_Cont == "Control" ~ "Control_PacBio",
    GRset_common$Dataset == "Array" & GRset_common$Case_Cont == "Case" ~ "Test Kabuki",
    GRset_common$Case_Cont == "Case" ~ "Test Kabuki",
    GRset_common$Case_Cont == "Control" ~ "Test Control"
  ),
  levels = c(
    "Train Kabuki", "Train Control", "Test Kabuki",
    "Kabuki_ONT_A", "Kabuki_ONT_WGS", "Kabuki_PacBio",
    "Test Control", "Control_ONT_A", "Control_ONT_WGS", "Control_PacBio",
    "CHARGE"
  )
)

GRset_common$dataset_id <- GRset_common$Dataset
GRset_common$dataset_id[GRset_common$dataset_id %in% c("PacBio", "ONT_A", "ONT_WGS", "Array")] <- "PMID41225292"
GRset_common$dataset_id <- factor(
  GRset_common$dataset_id,
  levels = c("GSE97362", "GSE116300", "GSE218186", "PMID41225292")
)
```

## PCA on Kabuki episignature CpGs


``` r
betas <- t(minfi::getBeta(GRset_common[model_cpgs, ]))

pca_res <- stats::prcomp(betas, scale. = TRUE)
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)

pca_df <- as.data.frame(pca_res$x) %>%
  dplyr::select(PC1, PC2) %>%
  tibble::rownames_to_column(var = "geo_accession") %>%
  dplyr::left_join(as.data.frame(SummarizedExperiment::colData(GRset_common)), by = "geo_accession")

col_pal <- c(
  "Train Control" = "lightblue",
  "Train Kabuki" = "lightcoral",
  "Test Control" = "dodgerblue",
  "Control_ONT_A" = "cornflowerblue",
  "Control_ONT_WGS" = "royalblue",
  "Control_PacBio" = "mediumblue",
  "Test Kabuki" = "firebrick",
  "Kabuki_ONT_A" = "orangered",
  "Kabuki_ONT_WGS" = "#7a2607",
  "Kabuki_PacBio" = "darkred",
  "CHARGE" = "forestgreen"
)

shape_pal <- c(
  "GSE97362" = 21,
  "GSE116300" = 22,
  "GSE218186" = 23,
  "PMID41225292" = 24
)

plt_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2)) +
  ggplot2::geom_point(ggplot2::aes(fill = Status, shape = dataset_id), size = 2.5, alpha = 0.75, color = "black") +
  ggplot2::labs(
    x = paste0("PC1 (", var_explained[1], "% var)"),
    y = paste0("PC2 (", var_explained[2], "% var)"),
    fill = "Status",
    shape = "Dataset"
  ) +
  ggplot2::scale_fill_manual(values = col_pal) +
  ggplot2::scale_shape_manual(values = shape_pal) +
  plot_custom_theme()
```

## SVM scores across datasets


``` r
predictions_df <- predict_model_on_grset(svm_mod, GRset_common) %>%
  dplyr::left_join(as.data.frame(SummarizedExperiment::colData(GRset_common)), by = "geo_accession")

plt_svm <- svm_plt(predictions_df, col_pal = col_pal, shape_pal = shape_pal)

openxlsx::write.xlsx(
  predictions_df,
  file.path(path_save, "Tables", "SVM_predictions_Episignature.xlsx"),
  rowNames = FALSE
)
```

The score plot answers a classification question:

- confirmed Kabuki array samples should receive high case scores if the expected episignature is preserved;
- controls should remain low if specificity is maintained;
- CHARGE/CHD7 samples test whether the classifier confuses a related developmental disorder with Kabuki;
- ONT and PacBio samples test whether the array-derived Kabuki signal remains detectable in long-read methylation estimates.

In the manuscript analysis, the molecularly confirmed Kabuki samples and controls were classified according to their expected status, and variant samples previously interpreted as non-pathogenic produced negative results consistent with the published interpretation. The long-read Kabuki and control profiles were also classified as expected. These observations demonstrate transfer of a strong signal in this example, not universal cross-platform equivalence.

## Long-read deviation analysis

The final part of `Kabuki_analysis.R` asks a different question: even if classification is correct, how far are the methylation values from an array-based reference? The script:

1. selects high-confidence non-long-read samples: cases with `Case > 0.75` and controls with `Case < 0.25`;
2. randomly uses 80% of those samples to build class-specific array reference profiles;
3. calculates the mean beta-value and beta-value standard deviation for each CpG within cases and controls;
4. calculates each remaining sample's deviation from the appropriate class reference;
5. compares those deviations across array and long-read platforms.

These 0.75/0.25 thresholds are used here to construct a clean reference set for this exploratory comparison. They are not proposed as general diagnostic thresholds.


``` r
set.seed(123)

samples_keep <- predictions_df %>%
  dplyr::filter(
    dataset_id != "PMID41225292",
    Status != "CHARGE",
    (Case > 0.75 & Case_Cont == "Case") | (Case < 0.25 & Case_Cont == "Control")
  ) %>%
  dplyr::group_by(Status, Case_Cont) %>%
  dplyr::slice_sample(prop = 0.8) %>%
  dplyr::pull(geo_accession)

GRset_reference <- GRset_common[, colnames(GRset_common) %in% samples_keep]

summary_betas <- as.data.frame(minfi::getBeta(GRset_reference)) %>%
  tibble::rownames_to_column(var = "id") %>%
  tidyr::pivot_longer(-id, names_to = "geo_accession", values_to = "Beta") %>%
  dplyr::mutate(Beta = as.numeric(Beta)) %>%
  dplyr::left_join(as.data.frame(SummarizedExperiment::colData(GRset_reference)), by = "geo_accession") %>%
  dplyr::group_by(id, Case_Cont) %>%
  dplyr::reframe(
    mean_beta = mean(Beta, na.rm = TRUE),
    sd_beta = stats::sd(Beta, na.rm = TRUE),
    n = dplyr::n()
  )
```

Then compute sample-level deviations:


``` r
GRset_valid <- GRset_common[
  ,
  !colnames(GRset_common) %in% samples_keep &
    GRset_common$Status != "CHARGE"
]

diff_betas <- as.data.frame(minfi::getBeta(GRset_valid)) %>%
  tibble::rownames_to_column("id") %>%
  tidyr::pivot_longer(-id, names_to = "geo_accession", values_to = "Beta") %>%
  dplyr::mutate(Beta = as.numeric(Beta)) %>%
  dplyr::left_join(as.data.frame(SummarizedExperiment::colData(GRset_valid)), by = "geo_accession") %>%
  dplyr::left_join(summary_betas, by = c("id", "Case_Cont")) %>%
  dplyr::mutate(
    beta_diff = Beta - mean_beta,
    z_score = beta_diff / sd_beta
  )

sample_deviation <- diff_betas %>%
  dplyr::group_by(geo_accession, Case_Cont, dataset_id, Status) %>%
  dplyr::reframe(
    mean_abs_diff = mean(abs(beta_diff), na.rm = TRUE),
    mean_abs_z = mean(abs(z_score), na.rm = TRUE),
    mean_z = mean(z_score, na.rm = TRUE),
    sd_z = stats::sd(z_score, na.rm = TRUE),
    n_cpg = dplyr::n()
  ) %>%
  dplyr::left_join(predictions_df %>% dplyr::select(geo_accession, Case), by = "geo_accession")
```

This deviation analysis is exploratory. In the manuscript example, long-read profiles show larger deviations from the array-derived reference than independent array samples, while matched array profiles remain closer to the expected array range. This supports a platform-mediated shift in the methylation measurements even though the Kabuki disease signal is still detectable [@hildonen2026longread].

## What Kabuki teaches

The Kabuki workflow illustrates several validation questions that should be kept separate:

- **Does the fixed CpG panel distinguish Kabuki from controls in independent array data?**
- **Does it remain specific when challenged with CHARGE/CHD7 comparators?**
- **Can the disease-associated signal still be detected in ONT and PacBio methylation profiles?**
- **Are the underlying methylation values numerically comparable across platforms?**

For this Kabuki example, the answer to the classification question is encouraging, while the PCA, clustering, and deviation analyses still show systematic technology-related differences. The result is therefore best viewed as a proof of principle for a strong established episignature, not as evidence that every array-derived episignature can be transferred to long-read data without additional validation.



