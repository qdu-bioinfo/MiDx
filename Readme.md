# BioMetaPipe

## Content

- [Introduction](#introduction)
- [Package requirement](#package-requirement)
- [Installation environment](#installation-environment)
- [Data Description](#build-a-phylogenetic-tree)
  - [a. 16S rRNA sequence data](#a-16s-amplicon)
  - [b. WGS Metagenome data](#b-wgs-metagenome)
- [Hyper-parameter adjustment](#model-training-and-testing)
- [Five-fold cross validation](#run-the-example-with-one-click)
- [Leave-one-out cross validation](#run-the-example-with-one-click)
- [Independent Verification](#supplementary)
## Introduction
BioMetaPipe, a five-stage computational framework designed for microbial data analysis, covers modules such as data preprocessing, batch effect correction, feature screening and standardization modeling, and model validation.
## Package requirement
```
environment.yml
```
## Installation environment
```
1.https://github.com/qdu-bioinfo/MiDx.git
2.cd MiDx
3.conda env create -f environment.yml -n BioMetaPipe
4.conda activate BioMetaPipe
```
## Data Description
### a. 16S rRNA sequence data
For all 16S rRNA sequencing data, we used Parallel-Meta Suite  (PMS, v.3.73) based on the Greengenes2 database, with a 99% similarity threshold, and applied copy number correction to obtain taxonomic and functional profiles.

| SampleID | OTU_RS.GCF.016803185.1.NZ.WBJP01000152.1 | OTU_MJ006.2.barcode60.umi189791bins.ubs.4 | ...  |
| -------- | ------- | ------- |:----:|
| sample1  | 0.001   | 0.002   | ...  | 
| sample2  | 0       | 0.003   | ...  | 
| sample3  | 0.005   | 0       | ...  | 
### b. WGS Metagenome data
For all metagenomes, low-quality reads were removed using the Fastp  tool. The pre-filtered reads were then aligned to the human genome (hg19) using Bowtie2 for host removal. In the profiling step, taxonomic profiles were generated using MetaPhlAn v.4.0.6 (https://github.com/biobakery/MetaPhlAn)  with the reference database mpa_vOct22_CHOCOPhlAnSGB_202212. Functional profiling was conducted using HUMAnN v.3.7 (https://github.com/biobakery/humann) with reference database UniRef90s (version 201901b) and default parameters. The resulting abundance profiles are further summarized into other gene groups, such as KOs (KEGG Orthologs), by using the command “human_regroup_table”. 

| SampleID | s__Clostridia_bacterium | s__Rothia_mucilaginosa | ...  |
| -------- | ----------------------- | ---------------------- | :--: |
| sample1  | 0.001                   | 0.002                  |      |
| sample2  | 0                       | 0.003                  |      |
| sample3  | 0.005                   | 0                      |      |


## Hyper-parameter adjustment
The hyper-parameters were toned by BayesianOptimization of scikit-learn.
```
python ./src/Bayesian_Optimization.py
```
## Five-fold cross validation
```
python ./src/CV.py
```
## Leave-one-out cross validation
```
python ./src/LODO.py
```
Below is a description of all available parameters:

| Option                  | Type      | Required | Description                                                                                             |
|-------------------------|-----------|----------|---------------------------------------------------------------------------------------------------------|
| `--raw-data`            | string    | Yes      | Path to the raw feature table CSV (samples × features).                                                 |
| `--meta-data`           | string    | Yes      | Path to the sample metadata CSV (index = sample IDs).                                                    |
| `--feature-type`        | string    | Yes      | Data type tag (e.g., `16s`, `wgs`); selects feature input modality.                                      |
| `--groups`              | list      | Yes      | List of group labels to include (e.g., `CTR`, `CRC`, `ADA`).                                             |
| `--mapping`             | list      | Yes      | Mapping of each group label to integer codes (format `LABEL=INTEGER`).                                   |
| `--class-level`         | string    | No       | Classification level tag (e.g., `t_sgb`, `otu`); default is `t_sgb`.                                      |
| `--filter-mode`         | string    | No       | Prevalence filtering metric (e.g., `pass`, `abundance`); default is `pass`.                              |
| `--filter-threshold`    | float     | No       | Threshold for prevalence filtering; default is `0.0`.                                                     |
| `--norm-mode`           | string    | No       | Normalization method (e.g., `pass`, `std`, `log`); default is `pass`.                                     |
| `--correction`          | string    | No       | Batch-correction method (`MMUPHin`, `Combat`, `none` or custom tag); default is `MMUPHin`.               |
| `--feature-method`      | string    | No       | Feature statistics or selection method (e.g., `lefse`, `wilcoxon`); default is `lefse`.                   |
| `--output-path`         | string    | Yes      | Directory path where all outputs will be saved.                                                          |
| `--no-correction`       | flag      | No       | Skip the batch-correction stage.                                                                         |
| `--no-filter`           | flag      | No       | Skip filtering and normalization stages.                                                                  |
| `--no-norm`             | flag      | No       | Skip normalization stage.                                                                                |
| `--no-feature-selection`| flag      | No       | Skip the feature selection stage.                                                                        |
| `--no-cv`               | flag      | No       | Skip cross-validation modelling stage.                                                                   |
| `--no-lodo`             | flag      | No       | Skip the leave-one-dataset-out evaluation stage.                                                         | 

## Independent Verification
```
python ./src/Predict_independent.py
```
