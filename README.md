# SNU-Seq Data Processing and Analysis

## These scripts are used in the following studies:
1. [Gerlevik, U., Lorenz, P., Lamstaes, A., Fischl, H., Xi, S., Saukko-Paavola, A., Murray, S., Brown, T., Welch, A., George, C., Angel, A., Furger, A., Mellor, J. (2025). "Single Nucleotide Resolution 4sU Sequencing (SNU-Seq) reveals the transcriptional responsiveness of an epigenetically primed human genome". *bioRxiv*.](https://www.biorxiv.org/content/10.1101/2021.07.14.452379v2.full)
2. Gerlevik, U. (2025). "Understanding the Capabilities of SNU-Seq and Applying it to Investigate the Relationship between Transcription, Chromatin and Metabolism in mIDH Glioma Cells". *Ludwig Maximilian University of Munich*.

## Scripts in order

__1. "scripts/HEK293/1_genome_and_annotations"__
  - Build [STAR](https://github.com/alexdobin/STAR) genome with [hg38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) for mapping
  - Find genomic-A stretches and merge with [ENCODE blacklist](https://github.com/Boyle-Lab/Blacklist)
  - Prepare [GENCODE v46](https://www.gencodegenes.org/human/release_46.html) annotations for metagenes and counting transcription start site (TSS), gene body, gene end and readthrough regions
  - Prepare [housekeeping genes from RSeQC](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38.HouseKeepingGenes.bed.gz/download) for counting in quality control steps

__2. "scripts/HEK293/2_preprocessing_and_qualityControl"__
  - Read quality control via [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  - Adapter and quality trimming using [bbduk.sh](https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk.sh)
  - Map the reads to the built STAR genome
  - Report quality control via [MultiQC](https://github.com/MultiQC/MultiQC)
  - Principal component analysis via [deepTools](https://github.com/deeptools/deepTools/)

__3. "scripts/HEK293/3_filtering_normalisation_and_merging"__
  - Filter blacklist, genomic A stretches and top/bottom 0.5% signal as "outliers" skewing the data using [bedtools](https://github.com/arq5x/bedtools2) and awk
  - Count the spike-in reads via [featureCounts](https://doi.org/10.1093/bioinformatics/btt656)
  - Estimate size factors from spike-ins using [DESeq2](https://doi.org/10.1186/s13059-014-0550-8)
  - Merge replicates of the same libraries to enhance the signal via bedtools and awk

__4. "scripts/HEK293/4_subtacting_noPAP_from_bPAP"__
  - Count the reads at the 3' end region to normalise between +bPAP and no bPAP libraries using bedtools since they have different library complexity with an expectation of a fairly similar profile of polyadenylated transcripts at the 3' end
  - Estimate size factors from 3' end regions using DESeq2 and apply them via bedtools and awk
  - Subtract no bPAP from +bPAP at single nucleotide resolution using deepTools
  - Remove negatives emerged because of the unmatched regions between +bPAP and no bPAP, and quantify them using awk, [ggplot2](https://github.com/tidyverse/ggplot2) and [ggpubr](https://github.com/kassambara/ggpubr)
  - Generate bedgraphs to view on [Integrative Genomics Viewer (IGV)](https://igv.org) with negative values for the reverse strand using awk

__5. "scripts/HEK293/5_metagene_analysis"__
  - Compute and plot metagenes using deepTools

__6. "scripts/HEK293/6_spliceJunction_analysis"__
  - Calculate splicing efficiency using [SPLICE-q](https://github.com/vrmelo/SPLICE-q)
  - Calculate mean splicing efficiency per sample and visualise the distributions using ggplot2 and ggpubr

__7. "scripts/HEK293/7_synthesis_decay_pausing_rate_analysis"__
  - Count the reads at TSS, gene body, gene end and readthrough regions using bedtools map by summing the signal in the processed bedgraph files
  - Prepare a counts table involving all region counts and all samples and normalise the summed signal dividing by the width of the regions
  - Calculate synthesis & decay rates, pausing index and termination efficiency. k-means clustering the genes based on the synthesis rate and visualise them using ggpubr

__8. "scripts/HEK293/8_TSS_enrichment_in_bPAP_over_noPAP"__
  - Calculate the median fold change of the signal in +bPAP over no bPAP at the TSS region

__9. "scripts/HEK293/9_comparison_with_TTseq_PROseq_NETseq"__
  - Put [Phil's TT-Seq GSM5452296](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5452296) and publicly available [TT-Seq GSM4730176](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4730176), [PRO-Seq GSM4730174](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4730174) and [mNET-Seq GSM7990390](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7990390) data in HEK293 cells to the same scale
  - Compute and plot metagenes of SNU-Seq, TT-Seq, PRO-Seq and mNET-Seq using deepTools
  - Generate negative reverse strand bedgraph to visualise Phil's TT-Seq data on IGV

__10. "scripts/HEP3B/1_prepare_ATACseqPeaks_and_FANTOM5_forMetagenes"__
  - Prepare [Anna's ATAC-Seq GSE172053](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172053) peaks and [FANTOM5](https://fantom.gsc.riken.jp/5/) enhancers data in a similar way to the genome annotation preparation in HEK293 cells to get clean/non-intersecting regions for a reliable metagene representation

__11. "scripts/HEP3B/2_prepare_SNUseq_ChIPseq_data"__
  - Prepare [Anna's SNU-Seq and ChIP-Seq of H3K27ac and H3K4me3 data from GSE172053](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172053) by file type conversions, concatenation and log2 transformation
  - Calculate H3K27ac to H3K4me3 ratio using deepTools and summarise the resulting ratios

__12. "scripts/HEP3B/3_signal_presence_in_ATACpeak_and_FANTOM5_regions"__
  - Determine the SNU-seq, H3K27ac and H3K4me3 signal and quantify the distribution of them on the ATAC-Seq peaks and FANTOM5 regions

__13. "scripts/HEP3B/4_metagene_analysis"__
  - Sort the ATAC-Seq peaks and FANTOM5 annotations accordingly, and compute and plot metagenes using deepTools


## Dependencies & Environments

This pipeline utilizes multiple Conda environments to manage dependencies for different stages of the analysis. All environment configuration files are located in the [`envs/`](https://github.com/ugerlevik/SNU-seq/tree/main/envs) directory.

### Environment List

| Analysis Stage | Environment File | Key Tools |
| :--- | :--- | :--- |
| **DeepTools** | [`envs/deeptools_env.yml`](envs/deeptools_env.yml) | DeepTools, Python 3.8 |
| **QC aggregation** | [`envs/multiqc_env.yml`](envs/multiqc_env.yml) | MultiQC |
| **Genomic-A flagging** | [`envs/py27.yml`](envs/py27.yml) | Python 2.7 scripts |
| **Splicing analysis** | [`envs/spliceQ_env.yml`](envs/spliceQ_env.yml) | SPLICE-q |

### Installation

To replicate a specific environment, use the following command structure:

```bash
# Example: Creating the DeepTools environment
mamba env create -f envs/deeptools_env.yml

# Activate the environment
mamba activate deeptools_env
```

### R Dependencies
R packages are included within the respective conda YAML files where possible. For most of the R and package versions used in the scripts, please refer to [`envs/R_versions.txt`](https://github.com/ugerlevik/SNU-seq/tree/main/envs/R_versions.txt).

## References
1. [manschmi/MexNab_3seq](https://github.com/manschmi/MexNab_3seq)
2. [manschmi/MS_Metagene_Tools](https://github.com/manschmi/MS_Metagene_Tools)
3. [nf-core/rnaseq](https://github.com/nf-core/rnaseq)
4. See the used R packages and other tools in the scripts.
