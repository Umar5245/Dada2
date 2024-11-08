# Dada2
DADA2 Pipeline for Microbial Community Analysis

This document describes the implementation of the DADA2 pipeline for processing and analyzing high-throughput sequencing data. The workflow includes quality control, error learning, sequence table generation, taxonomic assignment, and visualization of microbial diversity.

Requirements

Before running the pipeline, ensure the following R packages are installed: dada2, BiocManager, phyloseq, Biostrings, and ggplot2.

To install required packages, use the following commands:

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
s\
Setup and Input Data

Set the path variable to the directory containing your sequencing files.
Ensure files follow the naming convention: SAMPLENAME_R1_001.fastq for forward reads and SAMPLENAME_R2_001.fastq for reverse reads.

Example:
path <- "/path/to/MiSeq_SOP"
list.files(path)

Quality Control

Visualize the quality of forward and reverse reads to guide trimming parameters.

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
Filtering and Trimming

Filter low-quality reads and save results to a filtered/ subdirectory.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
Error Model and Denoising

Learn error rates and apply the DADA2 algorithm to identify Amplicon Sequence Variants (ASVs).

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
Merging and Chimera Removal

Merge paired-end reads and remove chimeric sequences.

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
Taxonomic Assignment

Assign taxonomy using the Silva reference database.

taxa <- assignTaxonomy(seqtab.nochim, "path/to/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
Phyloseq Object Creation

Create a phyloseq object for downstream analysis.

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa))
Data Visualization

Generate richness plots and Bray-Curtis NMDS ordinations.


plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
Create a bar plot for the top 20 taxa.


plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")

Export Data

Save processed feature tables and metadata as CSV files.

write.csv(feature_table, "feature_table.csv", row.names=TRUE)
write.csv(metadata, "metadata.csv", row.names=TRUE)
File Structure

Input Files:

Fastq files (_R1_001.fastq, _R2_001.fastq)
Reference taxonomy database (silva_nr_v132_train_set.fa.gz)
Output Files:

taxonomy_assignments.csv: Taxonomic assignments of ASVs.
feature_table.csv: Feature table with ASV abundances.
metadata.csv: Metadata associated with the samples.