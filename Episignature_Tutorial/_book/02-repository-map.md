# Repository Map {#repo-map}

This chapter maps the analysis scripts used throughout the tutorial. It is intended as a simple reference: for each dataset, it shows which script prepares the data, which outputs are created, and where those outputs are used later.

## Top-level helper

`scripts/functions_TutorialEpis.R` contains shared functions:

- `read_idat_minfi()` reads IDATs, runs sample/probe detection P-value filters, predicts sex, normalizes, drops SNP and cross-reactive probes, and optionally estimates immune cell composition.
- `detect_pca_outliers()` identifies PC1/PC2 outliers by Mahalanobis distance within `Case_Cont` groups.
- model wrappers such as `svm_training_f()`, `elastic_net_training_f()`, `ranger_training_f()`, `knn_training_f()`, and others train binary classifiers through `caret`.
- plotting helpers create SVM score plots, heatmaps, and palettes.

The tutorial helper at `scripts/bookdown_episignatures_tutorial/R/functions_tutorial_epis.R` follows the same general analysis logic but makes paths explicit and adds input checks that are useful for a public, step-by-step workflow. These safeguards also address several portability and first-run issues documented in the maintainer notes.

## Williams scripts

### `scripts/Williams_analysis.R`

This is the main Williams episignature discovery workflow.

Its sections are:

1. configure paths for `GSE119778`;
2. retrieve phenotype metadata from GEO;
3. read IDATs and preprocess methylation arrays;
4. estimate immune cell proportions;
5. run PCA for batch and outlier detection;
6. define discovery and validation samples;
7. run limma on discovery samples;
8. rank DMPs and select episignature CpGs;
9. train several classifiers on beta-values and M-values;
10. score discovery, internal validation, and external `GSE66552` samples;
11. export plots and prediction tables.

Important outputs include:

- `Williams/GSE119778/GRset_GSE119778.rds`
- `Williams/GSE119778/Tables/DMPs_limma_GSE119778.xlsx`
- `Williams/GSE119778/Tables/EpisignatureCpGs_GSE119778.xlsx`
- `Williams/GSE119778/_classification_models_GSE119778.rds`
- `Williams/GSE119778/Tables/Classification_predictions_GSE119778.xlsx`

### `scripts/Williams/GSE66552/GSE66552_analysis.R`

This script prepares the independent Williams validation/comparator dataset `GSE66552` [@strong2015williams]. It:

1. parses GEO metadata for typical controls, Williams syndrome, and 7q11.23 duplication samples;
2. preprocesses the IDATs;
3. runs PCA to inspect technical and diagnostic structure;
4. subsets the object to the Williams CpGs discovered only from `GSE119778`;
5. saves `GRset_epi_GSE66552.rds` for external scoring by the main Williams script.

The 7q11.23 duplication samples are comparators rather than Williams cases and are useful for checking whether the learned methylation pattern is specific to the deletion-associated Williams profile.

## Kabuki scripts

### `scripts/Kabuki/GSE97362/Kabuki_GSE97362_analysis.R`

This script prepares the array cohort used to train the Kabuki classifier. Unlike the Williams workflow, it does not discover a new CpG panel. It reads an established Kabuki episignature panel stored at:

```text
scripts/Kabuki/PMID41225292/episignature_PMID41225292.xlsx
```

The script:

1. parses GEO metadata including sex, age, diagnosis, and sample type;
2. preprocesses IDATs with `read_idat_minfi()`;
3. subsets the data to the published Kabuki episignature CpGs;
4. identifies the confirmed discovery/training subset defined by the original Kabuki study;
5. trains a linear SVM on that fixed CpG panel;
6. scores the remaining samples and exports the prediction results.

The underlying Kabuki methylation cohort was described by Butcher et al. [@butcher2017chargekabuki].

### `scripts/Kabuki/GSE116300/Kabuki_GSE116300_analysis.R`

This prepares the independent Kabuki cohort `GSE116300`, originally reported by Sobreira et al. [@sobreira2017kabuki]. The script uses `process_meth = "quantile"` rather than functional normalization and compares IDAT-derived beta-values with the processed GEO values by PCA. This comparison helps the reader see whether the two preprocessing sources produce similar large-scale structure; it does not make the two preprocessing methods technically identical.

### `scripts/Kabuki/GSE218186/Kabuki_GSE218186_analysis.R`

This prepares the additional external Kabuki array dataset `GSE218186`. Reported sex is not available in the GEO metadata used by the script, so predicted sex is stored in the harmonized `GRset_epis` object. The dataset is used for external testing rather than for defining the Kabuki CpG panel.

### `scripts/Kabuki/PMID41225292/LR_analysis.R`

This script converts the methylation matrix provided with the long-read Kabuki study into a `GenomicRatioSet` using `makeGenomicRatioSetFromMatrix()`. It reads metadata rows such as `Platform`, `Affected`, and `Cohort`, creates the sample labels, and saves `GRset_epis_PMID41225292.rds` [@hildonen2026longread].

This object is the bridge between the array-derived Kabuki episignature and methylation estimates generated with ONT and PacBio, together with the matched array measurements included in that study.

### `scripts/Kabuki/Kabuki_analysis.R`

This is the combined Kabuki validation script. It:

1. loads the SVM trained in `GSE97362`;
2. loads the prepared `GRset_epis_*` objects;
3. identifies CpGs shared across the included datasets;
4. checks that the CpGs required by the trained model are available before scoring;
5. combines the compatible objects into `GRset_common`;
6. labels training samples, external array samples, CHARGE/CHD7 comparators, and ONT/PacBio samples;
7. runs PCA on the Kabuki episignature CpGs;
8. applies the trained SVM;
9. builds heatmap and summary figures;
10. compares long-read methylation values with high-confidence array reference profiles.

This script therefore addresses two questions separately: whether the Kabuki disease signal remains classifiable, and whether the numerical methylation values are directly comparable across technologies.

