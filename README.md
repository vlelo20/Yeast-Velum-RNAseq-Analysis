# **Transcriptomic Profiling of *Saccharomyces cerevisiae* During Velum Biofilm Formation**
Author: Vian Lelo
Date created: January 20th, 2026, Last Updated: February 8th, 2026

# **1.0 Introduction**
Saccharomyces cerevisiae flor strains form a structured floating biofilm, termed a velum, at the air–liquid interface during the biological aging of sherry-style wines. This aerobic lifestyle enables oxidative metabolism under conditions of low nutrient availability, elevated ethanol concentration, and oxidative stress, driving the accumulation of acetaldehyde, acetals, and other flavor-active compounds characteristic of fino and manzanilla wines (Moreno-García et al., 2017). Velum development progresses through distinct morphological stages—early, thin, and mature—and prior transcriptomic and proteomic work has demonstrated stage-specific remodeling of gene expression, including strong induction of the cell-adhesion gene FLO11 and broad activation of stress-response and respiratory programs during maturation (Mardanov et al., 2020; Moreno-García et al., 2017).

The present study analyzes bulk RNA-seq data from BioProject PRJNA592304, spanning three velum stages with ≥3 biological replicates per stage, with three objectives: (i) robustly quantify transcript abundance, (ii) identify differentially expressed genes (DEGs) across stage-to-stage contrasts, and (iii) perform functional enrichment analysis to connect transcriptomic changes with the biological processes underlying velum maturation.

Transcript quantification will be performed using Salmon in selective-alignment, decoy-aware mode against the S. cerevisiae Ensembl R64-1-1 reference. Selective alignment improves mapping sensitivity over lightweight k-mer hashing alone, while genomic decoys reduce spurious assignment of reads from unannotated loci—both are meaningful considerations in the compact yeast genome (Patro et al., 2017). Raw read quality will be assessed with FastQC and aggregated using MultiQC to detect per-sample outliers prior to downstream analysis (Andrews, 2010; Ewels et al., 2016). Differential expression will be tested with DESeq2, which employs negative-binomial generalized linear models with median-of-ratios normalization and empirical Bayes shrinkage—a framework well-suited to small replicate designs (Love et al., 2014). Results will be cross-validated using edgeR, which applies trimmed mean of M-values (TMM) normalization and quasi-likelihood testing, offering complementary modeling assumptions and serving as a sensitivity check (Robinson et al., 2010). Where appropriate, limma-voom will additionally be considered; by applying precision weights to log-counts per million, it bridges count-based and linear model frameworks and performs comparably to negative-binomial methods at moderate sample sizes (Law et al., 2014). Finally, over-representation analysis and/or gene-set enrichment analysis will be applied to ranked DEG lists to contextualize findings within established pathways of adhesion, respiration, and stress response relevant to velum biology.

Each tool carries trade-offs. Salmon does not produce genome-aligned BAM files, precluding visualization of novel splice junctions—a minor limitation in the well-annotated yeast genome but relevant to consider. DESeq2's shrinkage estimators stabilize variance at low counts but assume independence across genes; edgeR's quasi-likelihood pipeline relaxes distributional assumptions but can be sensitive to filtering choices; and limma-voom's linear model framework is less appropriate when count distributions are highly overdispersed. Together, these methods provide a robust, multi-framework approach to characterizing the transcriptional architecture of flor yeast velum development.

# **2.0 Methods**
## **2.1 Environment Setup**
To ensure a standardized and reproducible computational environment, all analyses were performed within a dedicated Conda environment (Python 3.10) containing sra-tools, fastqc, multiqc, fastp, and salmon. A systematic directory structure was established to segregate raw data, trimmed reads, reference files, and logs. Final transcript quantification results and quality control metrics were consolidated into categorized output directories to maintain data provenance throughout the pipeline.
```
conda create -n transcriptomics python=3.10 sra-tools pigz fastqc multiqc fastp salmon -c bioconda -c conda-forge
conda activate transcriptomics

#version checks
salmon --version   
fastqc --version
fasterq-dump --version

#directory structure
mkdir -p ~/binf6110/assignment2/{raw_fastq,trimmed_fastq,reference,results/{fastqc_raw,fastqc_trim,multiqc,salmon_quant},logs}
cd ~/binf6110/assignment2
<img width="925" height="392" alt="image" src="https://github.com/user-attachments/assets/f025560c-ca1f-4334-a11a-0c035e3b8b48" />

```

## **2.2 Data & Metadata**

### **2.2.1 Dataset**
```
```
### **2.2.2 Metadata File**
```
```
### **2.2.3 Data Download**
**Raw Sequencing Data Acquisition**
Raw sequencing reads for nine samples—representing early (IL20–IL22), thin (IL23–IL25, IL29), and mature (IL30–IL31) biofilm developmental stages—were retrieved from the NCBI Sequence Read Archive (SRA). Data acquisition was performed using the SRA Toolkit utility fasterq-dump (v3.2.1), utilizing a manifest-driven approach to ensure experimental reproducibility and metadata synchronization. To facilitate downstream analysis, the original SRR accessions were programmatically mapped and renamed to their corresponding biological sample identifiers. Initial sequence quality was assessed using FastQC (v0.12.1) to ensure data integrity; for a complete record of the implementation parameters, please refer to the download.sh script in the scripts file.

**Reference Genome and Annotation Retrieval**
The Saccharomyces cerevisiae R64-1-1 reference assembly and its corresponding functional annotations were obtained from the Ensembl database (Release 110). Specifically, the cDNA sequences, the top-level genomic DNA assembly, and the Gene Transfer Format (GTF) annotation files were downloaded using the wget utility. These genomic resources provided the necessary framework for subsequent transcript quantification and analysis. To ensure high-performance mapping, the reference transcriptome was indexed using Salmon (v1.10.3). 

```
# From Ensembl FTP mirrors (change release if needed)
wget -P reference/ https://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
wget -P reference/ https://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
wget -P reference/ https://ftp.ensembl.org/pub/release-110/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz
```

## **2.3 Quality Control**
```
```
### **2.3.1 FastQC & MultiQC**

Following the data acquisition phase, a comprehensive quality control (QC) assessment was performed to evaluate the integrity and sequencing metrics of the raw reads. This was executed by running FastQC (v0.12.1) on all compressed FASTQ files within the project directory, utilizing 16 computational threads to expedite the processing of the nine biofilm samples. To synthesize the individual sample reports into a single, interpretable diagnostic overview, MultiQC was employed to aggregate the results into a unified HTML report. This dual-stage QC pipeline allowed for the systematic identification of potential sequencing artifacts, such as adapter contamination or base-calling quality drops, prior to downstream transcriptomic analysis; for specific implementation details, please refer to the quality control section in the scripts directory.
```
cd ~/binf6110/assignment2
fastqc data/fastq/*.fastq.gz -o results/fastqc_raw -t 16
multiqc results/fastqc_raw -o results/multiqc
```
### **2.3.2 Optional Trimming**
```
```
## **2.4 Reference & Indexing**

To ensure high-performance mapping and reduce spurious alignments from genomic contamination, a "gentrome" index was constructed. This involved concatenating the transcriptome and the whole-genome sequences, with the latter serving as a decoy. Genome targets were identified and extracted into a decoy list to guide the mapping process. The final selective-alignment index was generated using Salmon (v1.10.3) with a sensitive k-mer length, providing the computational framework for subsequent transcript quantification. For the exact source URLs and retrieval parameters used in this step, please refer to the reference download section within the scripts directory.
```
cd reference

# make a 'decoys.txt' from the genome FASTA headers
zcat Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz \
  | grep '^>' | cut -d ' ' -f1 | sed 's/>//' > decoys.txt

# build the "gentrome" (transcriptome + genome concatenated; genome last)
cat Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
    Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz \
  > gentrome.fa.gz

# make the index (use threads as you like, e.g. 16 on your 22‑thread box)
salmon index -t gentrome.fa.gz -d decoys.txt -p 16 -i salmon_index
cd ..
```
## **2.5 Quantification**
```
```
## **2.6 Differential Expression Analysis**
```
```
### **2.6.1 PCA & Sample Distances**
```
```
### **2.6.2 Pairwise Contrasts**
```
```
### **2.6.3 LRT Analysis**
```
```
## **2.7 Functional Annotation & Enrichment**
```
```
### **2.7.1 ORA**
```
```
### **2.7.2 GSEA**
```
```

# **3.0 Results**
## **3.1 Overall Data Structure**
## **3.2 Differential Expression**
## **3.3 Functional Annotation**

# **4.0 Discussion**
## **4.1 Biological Interpretation**
## **4.2 Implications**
## **4.3 Limitations**
## **4.4 Future Directions**

# **5.0 Reproducibility**
## **5.1 Version Information**
## **5.2 Pipeline Execution**

# **6.0 Repository Structure**



# **7.0 How to Run**

# **8.0 References**
Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data. Babraham Institute. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). Voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology, 15, R29. https://doi.org/10.1186/gb-2014-15-2-r29
Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550. https://doi.org/10.1186/s13059-014-0550-8
Mardanov, A. V., Eldarov, M. A., Beletsky, A. V., Tanashchuk, T. N., Kishkovskaya, S. A., & Ravin, N. V. (2020). Transcriptome profile of yeast strain used for biological wine aging revealed dynamic changes of gene expression in course of flor development. Frontiers in Microbiology, 11, 538. https://doi.org/10.3389/fmicb.2020.00538
Moreno-García, J., Mauricio, J. C., Moreno, J., & García-Martínez, T. (2017). Differential proteome analysis of a flor yeast strain under biofilm formation. International Journal of Molecular Sciences, 18(4), 720. https://doi.org/10.3390/ijms18040720
Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14(4), 417–419. https://doi.org/10.1038/nmeth.4197
Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: A Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616
