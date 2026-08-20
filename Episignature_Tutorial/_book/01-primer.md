# Primer for New Readers {#primer}

This chapter gives the minimum background needed to follow the workflow. It is written for readers who may be new to DNA methylation analysis, so the emphasis is on what each step does, why it is needed, and how it connects to episignature discovery or testing.

## Clinical epigenomics in rare disease

Rare diseases and neurodevelopmental disorders can remain difficult to diagnose even after extensive genomic testing. A patient may have no clearly causal variant, or testing may identify a variant of uncertain significance (VUS) whose functional effect is unclear. DNA methylation profiling provides a complementary readout because it can capture downstream molecular consequences of disease rather than the DNA sequence alone. In a growing number of rare disorders, reproducible disease-associated methylation profiles can be used as episignatures to support diagnosis and variant interpretation [@arefeshghi2020episignatures; @kerkhof2024reporting].

This tutorial was prepared by the **Bioinformatics Research Group in Genetic Epidemiology, ISGlobal** as a practical companion to the manuscript chapter on episignatures in rare disorders. It uses public datasets to demonstrate two complementary scenarios:

1. **Williams syndrome**: deriving a de novo episignature from methylation-array data, selecting CpGs, training classifiers, and testing external samples.
2. **Kabuki syndrome**: starting from an established CpG panel and evaluating whether array-trained models behave consistently across array and long-read methylation data.

## Glossary

**Episignature**: a reproducible, disease-associated DNA methylation profile. In practice, it is usually represented by a selected panel of CpG sites whose combined pattern separates affected individuals from appropriate reference samples. The diagnostic information comes from the multivariate pattern, not necessarily from any single CpG.

**Epivariant, epimutation, or rare differentially methylated region**: a localized methylation abnormality affecting a specific genomic region. This is different from a genome-wide episignature. Imprinting defects are a special case because they involve abnormal methylation at parent-of-origin-specific regulatory regions and may not present as a broad episignature.

**Classifier score**: a numerical output from a machine-learning model. In this tutorial, the main score is the model-derived probability assigned to the `Case` class. A high or low score can support interpretation, but it must be considered together with phenotype, genotype, sample quality, platform, and the reference data used to train the model.

## DNA methylation in one page

DNA methylation is an epigenetic modification in which a methyl group is added to cytosine, most commonly at **CpG sites**. Methylation patterns can be stably maintained through cell division while also changing with development and biological context. DNA methylation contributes to gene regulation, cell identity, X-chromosome inactivation, genomic imprinting, and repression of repetitive elements [@jones2012functions; @moore2013dna].

For the beta-values used in this tutorial, methylation is represented on a scale from 0 to 1:

- 0 means that little or no methylation is detected at that CpG.
- 1 means that methylation is detected in almost all measured molecules at that CpG.

Blood is a mixture of cell types, and each cell type has its own methylation profile. Most clinical episignature studies use peripheral blood because it is accessible and supported by large reference datasets. A negative blood-based result, however, does not prove that disease-relevant methylation changes are absent in another tissue [@kerkhof2024reporting].

Illumina methylation arrays measure hundreds of thousands of predefined CpG sites across the genome. The raw scanner files are **IDATs**, with one red-channel file and one green-channel file per sample. The R package `minfi` can read these files into an `RGChannelSet`, perform preprocessing and normalization, map probes to the genome, and produce a `GenomicRatioSet` (`GRset`) containing methylation measurements [@aryee2014minfi].

Two transformations are used throughout the scripts:

- **Beta-values**: values between 0 and 1 that are easy to interpret as relative methylation levels. In this tutorial they are used for visualization and for most classifier inputs.
- **M-values**: a logit transformation of beta-values. They are less intuitive to plot, but they have more suitable statistical properties for differential methylation modelling, particularly near beta-value extremes [@du2010comparison].

## What is an episignature?

An **episignature** is a reproducible disease-associated DNA methylation pattern that can function as a molecular biomarker. In this repository, the term refers to a selected panel of CpGs whose combined methylation profile can distinguish cases from controls and can be used in both unsupervised visualization and supervised classification [@arefeshghi2020episignatures].

There are two related but different tasks:

1. **Discovery** asks: which CpGs differ between cases and controls?
2. **Detection** asks: given a new sample, does its methylation profile look like the case group?

The scripts move from discovery to detection:

1. read and normalize IDAT files;
2. perform QC and remove problematic samples/probes;
3. estimate blood cell composition;
4. use PCA to inspect batch effects and outliers;
5. call differentially methylated positions with limma [@ritchie2015limma];
6. choose a compact, non-redundant CpG panel;
7. train a classifier, usually an SVM;
8. test held-out and independent samples.

## Cases, controls, and why matching matters

A classifier learns the differences present in its training data. If cases and controls also differ strongly in age, sex, blood cell composition, chip, laboratory, or array type, the model may learn those technical or biological differences instead of the disease-associated signal. The workflow therefore checks or adjusts for several potential confounders:

- sex prediction and reported-sex checks;
- removal of sex chromosome probes for discovery if necessary;
- estimated immune cell composition in the limma model;
- batch effect diagnosis using top variable CpGs, with variable importance for known confounders;
- held-out validation samples and external datasets.

These steps reduce risk, but they do not replace thoughtful cohort design. Feature selection and model fitting must also be performed using discovery/training samples only; using validation samples during those steps can produce overly optimistic performance estimates [@ambroise2002selection].

## Why blood cell composition is included

Blood is a mixture of cell types, and each cell type has a characteristic methylation profile. If cases and controls have different immune-cell proportions, those differences can appear as disease-associated methylation changes. Estimating cell composition and including appropriate cell fractions as covariates can reduce this source of confounding [@houseman2012cellmixture; @salas2018blood].

## Why PCA appears repeatedly

PCA is used for two distinct purposes:

1. **Technical QC**: before discovery, PCA is run on the most variable autosomal CpGs. The goal is to see whether slide, array position, or outlier samples dominate the data.
2. **Episignature visualization**: after CpG selection, PCA is rerun using only episignature CpGs. The goal is to see whether cases and controls separate in the selected feature space.

The first PCA is a quality-control view: it helps identify technical structure or unusual samples that need review. The second PCA is a signature view: it shows whether the selected CpGs separate the groups in an unsupervised projection. PCA is descriptive and does not by itself prove diagnostic performance.

## What an SVM score means

The `caret` SVM workflow used here returns class probabilities [@kuhn2008caret]. The main output column is called `Case` and ranges from 0 to 1:

- scores near 0 are more control-like;
- scores near 1 are more case-like;
- intermediate scores require more caution.

A threshold such as 0.5 is convenient for illustration, but it is not a universal clinical cutoff. Thresholds and score interpretation should be established with independent validation data and, where possible, disease-control cohorts.

## Discovery panel versus validation panel

The Williams workflow discovers its own CpG panel from `GSE119778` and then tests it in held-out samples and the independent `GSE66552` cohort. The Kabuki workflow is different: it starts from an established Kabuki CpG panel, stored in `scripts/Kabuki/PMID41225292/episignature_PMID41225292.xlsx`, and evaluates that fixed signal across independent array and long-read datasets [@arefeshghi2020episignatures; @hildonen2026longread]. This distinction is central to the tutorial:

- Williams teaches **how to discover and validate a new candidate signature**.
- Kabuki teaches **how to test an existing signature across cohorts and technologies**.


