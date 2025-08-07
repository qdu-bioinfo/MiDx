\name{NEWS}
\title{News for Package \pkg{tidybulk}}

\section{Changes in version 1.2.0, Bioconductor 3.12 Release}{
\itemize{
    \item Make gene filtering functionality `identify_abundance` explicit, a warning will be given if this has not been performed before the majority of workflow steps (e.g. `test_differential_abundance`).
    \item Add Automatic bibliography `get_bibliography`.
    \item Add DESeq2 and limma-voom to the methods for `test_differential_abundance` (method="DESeq2").
    \item Add prefix to test_differential_abundance for multi-methods analyses.
    \item Add other cell-type signature to `deconvolve_cellularity`.
    \item Add differential cellularity analyses `test_differential_cellularity`.
    \item Add gene descrption annotation utility `describe_transcript`.
    \item Add `nest` functionality for functional-coding transcriptomic analyses.
    \item Add gene overrepresentation functionality `test_gene_overrepresentation`.
    \item Add github website.
    \item Seep up data frame vadidation.
    \item Several bug fixes.
}}

\section{Changes in version 1.3.2, Bioconductor 3.13 Release}{
\itemize{
    \item Tidybulk now operates natively with SummarizedExperment data container, in a seamless way thanks to tidySummarisedExperiment 10.18129/B9.bioc.tidySummarizedExperiment
    \item Added robust edgeR as it outperforms many other methods as shown here doi.org/10.1093/nargab/lqab005
    \item Added test stratifiction cellularity, to easily calculate Kaplan-Meier curves
    \item Production of SummarizedExperiment from BAM or SAM files
    \item Added treat method to edgeR and voom differential transcription tests doi.org/10.1093/bioinformatics/btp053
    \item Added the method as_SummarizedExperiment
    \item Vastly improved test_gene_enrichment
    \item Added test_gene_rank, based on GSEA
    \item Several bug fixes.
}}

\section{Changes in version 1.5.5, Bioconductor 3.14 Release}{
\itemize{
    \item Added user-defined gene set for gene rank test
    \item Sped up aggregate_transcripts for large scale tibbles or SummarizedExperiment objects
    \item Allow passing additional arguments to DESeq2 method in test_differential_abundance
    \item Allow scale_abundance to run with a user-defined subset of genes (e.g. housekeeping genes)
    \item Add UMAP to reduce_dimensions()
    \item Several minor fixes, optimisations and documentation improvements
}}

\section{Changes in version 1.7.3, Bioconductor 3.15 Release}{
\itemize{
    \item Improve imputation and other features for sparse counts
    \item Cibersort deconvolution, check 0 counts
    \item Improve missing abundance with force scaling
    \item Other small fixes and messaging
}}

\section{Changes in version 1.7.4, Bioconductor 3.16 Dev}{
\itemize{
    \item Improved deconvolution robustness for SummarizedExperiment, edge cases
    \item Allow mapping of tidybulk_SAM_BAM to non-human genomes
    \item Adopt the vocabulary .feature, .sample, for conversion between SummarizedExperiment and tibble, similarly to tidySummarizedExperiment
    \item Deprecate .contrasts argument if favour of contrasts (with no dot)
    \item Make aggregate_duplicates more robust for tibble and SummarizedExperiment inputs
    \item Deprecate log_tranform argument for all methods for a more generic tranform argument that accepts arbitrary functions
}}

\section{Changes in version 1.9.2, Bioconductor 3.16 Dev}{
\itemize{
    \item Improve aggregate_duplicates for tibble and SummarizedExperiment
    \item Fix epic deconvolution when using DelayedMatrix
    \item Allow as_SummarizedExperiment with multiple columns identifiers for .sample and .feature
}}
