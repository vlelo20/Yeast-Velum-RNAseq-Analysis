# **Transcriptomic Profiling of *Saccharomyces cerevisiae* During Velum Biofilm Formation**
Author: Vian Lelo
Date created: March 1st, 2026, Last Updated: March 1st, 2026

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

Raw RNA-seq reads were obtained from the NCBI Sequence Read Archive under BioProject accession PRJNA592304 (Mardanov et al., 2020). The dataset comprises nine paired-end Illumina samples derived from Saccharomyces cerevisiae flor yeast at three stages of velum development: Early biofilm (IL20, IL21, IL22; SRR10551665–SRR10551663), Thin biofilm (IL23, IL24, IL25; SRR10551662–SRR10551660), and Mature biofilm (IL29, IL30, IL31; SRR10551659–SRR10551657), with three biological replicates per stage. Reads were downloaded using fasterq-dump from SRA-tools and quantified against the S. cerevisiae Ensembl R64-1-1 reference transcriptome (Saccharomyces_cerevisiae.R64-1-1.110) using Salmon in selective-alignment, decoy-aware mode.

### **2.2.2 Metadata File**
A sample metadata file (samplesheet.csv) was constructed manually to map each sample identifier to its corresponding SRA run accession and velum stage, containing three columns: sample (e.g., IL20), srr (e.g., SRR10551665), and stage (Early, Thin, or Mature). The stage variable was encoded as an ordered factor with levels Early < Thin < Mature to reflect the biological progression of velum development, and was used as the grouping variable in all downstream DESeq2 models. No additional physicochemical metadata (e.g., ethanol content, aldehyde concentrations) reported in the original publication were incorporated into the statistical models, as the analysis was restricted to transcriptomic stage comparisons using the RNA-seq data alone.

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

### **2.3.1 FastQC & MultiQC**

Following the data acquisition phase, a comprehensive quality control (QC) assessment was performed to evaluate the integrity and sequencing metrics of the raw reads. This was executed by running FastQC (v0.12.1) on all compressed FASTQ files within the project directory, utilizing 16 computational threads to expedite the processing of the nine biofilm samples. To synthesize the individual sample reports into a single, interpretable diagnostic overview, MultiQC was employed to aggregate the results into a unified HTML report. This dual-stage QC pipeline allowed for the systematic identification of potential sequencing artifacts, such as adapter contamination or base-calling quality drops, prior to downstream transcriptomic analysis; for specific implementation details, please refer to the scripts directory.
```
cd ~/binf6110/assignment2
fastqc data/fastq/*.fastq.gz -o results/fastqc_raw -t 16
multiqc results/fastqc_raw -o results/multiqc
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

Transcript abundance for the nine biofilm samples was quantified using Salmon (v1.10.3) in its mapping-based quantification mode. Single-end reads were processed by aligning against the previously constructed S. cerevisiae R64-1-1 selective-alignment ("gentrome") index. To account for library-specific artifacts, the Salmon automatic library type detection (-l A) was employed. Furthermore, quantification fidelity was enhanced through the inclusion of bias correction parameters—specifically --gcBias and --seqBias—which mitigate technical variations inherent in sequencing library preparation, while the --validateMappings flag was utilized to refine alignment sensitivity. Calculations were executed across 16 parallel threads per sample to optimize computational throughput. The resulting quantification output files (quant.sf) provided both raw counts and normalized abundance estimates (Transcripts Per Million; TPM) for each transcript, which served as the foundation for subsequent differential gene expression analysis across the early, thin, and mature biofilm developmental stages. please refer to the scripts directory.

## **2.6 Differential Expression Analysis**
Differential expression analysis was performed using DESeq2 (v1.42; Love et al., 2014), with transcript-level abundance estimates imported via tximport using the countsFromAbundance = "no" setting to retain raw inferred counts. A negative-binomial generalised linear model was fitted with velum stage (Early, Thin, Mature) as the sole explanatory variable. Genes with fewer than 10 counts in fewer than three samples were removed prior to modelling. Size factors were estimated using DESeq2's median-of-ratios method, and dispersion was estimated with empirical Bayes shrinkage. Log2 fold changes were further shrunk using the adaptive shrinkage estimator (ashr; Stephens, 2017) to reduce noise for low-count genes. A Benjamini-Hochberg false discovery rate (FDR) threshold of 0.05 and an absolute log2 fold-change threshold of 1 were applied to define differentially expressed genes (DEGs).

### **2.6.1 PCA & Sample Distances**
To assess overall transcriptomic structure and sample reproducibility, count data were transformed using the variance-stabilising transformation (VST) implemented in DESeq2 with blind = FALSE. Principal component analysis (PCA) was performed on the VST-normalised matrix and visualised using ggplot2, with sample labels added via ggrepel to facilitate identification of potential outliers. Sample-to-sample Euclidean distances were computed on the transposed VST matrix and displayed as a hierarchically clustered heatmap using pheatmap, with samples annotated by velum stage. These diagnostics were used to confirm that biological replicates clustered by stage prior to statistical testing.

### **2.6.2 Pairwise Contrasts**
Three pairwise contrasts were extracted from the fitted DESeq2 model to characterise transcriptional changes at successive stages of velum development: Thin vs Early, Mature vs Early, and Mature vs Thin. For each contrast, results were obtained using the results() function with an independent filtering alpha of 0.05. Volcano plots were generated for all three contrasts, displaying shrunk log2 fold changes against -log10(adjusted p-value), with the top ten up- and down-regulated genes labelled by gene identifier. A combined three-panel figure was produced using the patchwork package to facilitate direct visual comparison across contrasts. Full results tables were exported as CSV files for downstream interpretation.

### **2.6.3 Likelihood Ratio Test**
To identify genes exhibiting any significant variation in expression across the three velum stages, irrespective of the direction or stage of change, a likelihood ratio test (LRT) was performed by comparing the full model (~stage) against a reduced intercept-only model (~1). This approach detects genes with a stage-dependent expression profile without being restricted to a single pairwise comparison, and is particularly suited to ordered developmental time-series such as the Early-Thin-Mature progression. Genes were considered significant at FDR < 0.05, and the resulting gene set was used to generate a clustered heatmap of the top 50 stage-varying genes to visualise temporal expression patterns.

## **2.7 Functional Annotation & Enrichment**
Functional enrichment analyses were performed to contextualise DEGs within known biological processes relevant to velum formation, including cell adhesion, oxidative metabolism, and stress response. Gene identifiers used throughout were Saccharomyces Genome Database (SGD) open reading frame (ORF) designations, compatible with the org.Sc.sgd.db annotation package (Bioconductor).

### **2.7.1 Over-Representation Analysis (ORA)**
Over-representation analysis (ORA) was conducted using the enrichGO() function from clusterProfiler (Wu et al., 2021), testing for enrichment of Gene Ontology Biological Process (GO:BP) terms among DEGs from the Mature vs Early contrast. The background gene universe was defined as all genes passing the pre-filtering step. P-values were adjusted using the Benjamini-Hochberg method, and terms with adjusted p < 0.05 and q < 0.2 were considered significantly enriched. Results were visualised as a dot plot showing the top 20 enriched terms, with dot size proportional to gene-set size and colour indicating adjusted p-value.

### **2.7.2 Gene Set Enrichment Analysis (GSEA)**
Gene set enrichment analysis (GSEA) was performed using the gseGO() function from clusterProfiler to capture coordinated shifts in GO:BP terms without imposing a binary significance threshold. Genes were pre-ranked by the DESeq2 Wald statistic from the Mature vs Early contrast, which integrates both fold-change magnitude and statistical confidence. The analysis was run with 1,000 permutations and a minimum gene-set size of 10. Normalised enrichment scores (NES) and FDR-adjusted p-values were reported, and the top positively and negatively enriched pathways were visualised using enrichplot. GSEA complements ORA by detecting pathways enriched among moderately but consistently regulated genes that may fall below a strict fold-change cut-off.

# **3.0 Results**

## 3.1 Transcriptome-Wide Structure and Sample Quality
To assess global transcriptomic structure and confirm replicate reproducibility prior to differential expression testing, variance-stabilised (VST) count data were subjected to principal component analysis (PCA) and pairwise Euclidean distance profiling. PCA revealed clear separation of all three velum stages along the first two principal components, which together explained 95% of total variance (PC1: 75%, PC2: 20%; Figure 1). Early biofilm replicates (IL20–IL22) clustered tightly in the lower-left quadrant, Thin biofilm samples (IL23–IL25) in the upper-left quadrant, and Mature biofilm samples (IL29–IL31) in the right-hand region of the plot. The clean stage-wise separation, combined with tight within-group clustering, indicated high biological reproducibility and confirmed that stage is the dominant source of transcriptional variance in this dataset.
Sample-to-sample Euclidean distance analysis corroborated these findings (Figure 2). Hierarchical clustering of the distance matrix grouped samples by stage, with each group forming a clearly demarcated block of low pairwise distances along the diagonal. Between-stage distances were consistently larger than within-stage distances, and no outlier samples were identified. Together, Figures 1 and 2 provide strong evidence that the dataset is of sufficient quality for downstream differential expression analysis.

<img width="1800" height="1500" alt="PCA_stage" src="https://github.com/user-attachments/assets/daa2a05e-df54-4c12-a712-ddd0ad6c61a9" />
Figure 1. Principal Component Analysis of VST-normalised RNA-seq counts. Each point represents one biological replicate, coloured by velum stage (Early = blue, Thin = orange, Mature = green). PC1 (75%) and PC2 (20%) together capture 95% of transcriptome-wide variance. Sample labels were added using ggrepel.

<img width="2100" height="1800" alt="sample_distance_heatmap" src="https://github.com/user-attachments/assets/1dfd05e0-3678-43bf-982e-3699c695b9ff" />
Figure 2. Sample-to-Sample Euclidean Distance Heatmap. Distances were computed on the transposed VST matrix and visualised with hierarchical clustering (complete linkage). Darker shading indicates smaller distance (greater similarity). Stage annotation bars are shown above columns. Replicates cluster tightly by stage with no outliers detected.


## 3.2 Differential Expression Across Velum Development Stages
Differential expression analysis was performed using DESeq2 with ashr log2 fold-change (LFC) shrinkage, applying thresholds of |LFC| > 1 and FDR < 0.05. Across all three pairwise contrasts, large numbers of differentially expressed genes (DEGs) were identified, consistent with the extensive transcriptional remodelling reported for flor yeast velum formation (Mardanov et al., 2020). The Thin vs Early contrast yielded 840 DEGs (410 down-regulated, 430 up-regulated), reflecting the initial transcriptional shift upon biofilm initiation. The Mature vs Early contrast produced the largest number of DEGs (2,116 total; 961 down, 1,155 up), indicating a cumulative and progressive transcriptional divergence over velum development. The Mature vs Thin contrast identified 1,526 DEGs (649 down, 877 up), demonstrating that substantial transcriptional change continues between the thin and mature biofilm stages.

Volcano plots for all three contrasts are shown in Figure 3. Notably, several genes appeared as top-ranked DEGs across multiple contrasts. YGR088W, YGR087C, and YHR094C were among the most significantly down-regulated genes in both the Thin vs Early and Mature vs Early comparisons, suggesting suppression of these loci early in velum development that is maintained through maturation. Conversely, YNR073C, YNR071C, and YNR072W were consistently up-regulated in the Mature vs Early and Mature vs Thin contrasts, pointing to genes with specific roles in mature biofilm physiology. The Mature vs Early contrast additionally showed strong induction of YJR152W and YGL055W, genes associated with oxidative and respiratory metabolism, consistent with the known shift toward aerobic respiration during velum maturation.

<img width="4800" height="1650" alt="volcano_all_contrasts" src="https://github.com/user-attachments/assets/5fcde508-5578-4d2a-81a7-54618b13cd62" />
Figure 3. Volcano plots for all three pairwise DESeq2 contrasts. Panels show Thin vs Early (left), Mature vs Early (centre), and Mature vs Thin (right). Red points indicate significantly up-regulated genes, blue points down-regulated genes, and grey points non-significant genes (|LFC| > 1, FDR < 0.05). The top 10 up- and down-regulated genes by adjusted p-value are labelled in each panel. Dashed vertical lines mark |LFC| = 1; dashed horizontal line marks FDR = 0.05.

To visualise expression patterns among the most significantly changing genes, a hierarchically clustered heatmap of the top 50 DEGs from the Mature vs Early contrast is shown in Figure 4. Two principal expression clusters were evident. The first cluster comprised genes with high expression in Early and Thin biofilm stages but reduced expression in Mature samples, including YGR088W, YHR094C, and YGR087C. The second cluster contained genes progressively induced during velum maturation, most strongly expressed in the Mature stage, including members of the YNR071C–YNR073C locus. Thin biofilm replicates occupied an intermediate transcriptional position, co-clustering with Mature samples in the column dendrogram rather than with Early samples, consistent with the PCA result and indicating that the largest transcriptional shift occurs between the Early and Thin stages.

<img width="2400" height="2700" alt="heatmap_top50_mature_vs_early" src="https://github.com/user-attachments/assets/8ff27d00-0603-4a1e-a3a9-5d4b455dc68c" />
Figure 4. Heatmap of the top 50 DEGs from the Mature vs Early contrast (row-centred VST counts). Rows represent genes ranked by adjusted p-value; columns represent individual samples annotated by stage. Colour scale ranges from deep blue (low relative expression) to deep red (high relative expression). Hierarchical clustering used Ward.D2 linkage.

3.3 Functional Enrichment Analysis
To identify biological processes enriched among the DEGs from the Mature vs Early contrast, over-representation analysis (ORA) was performed against Gene Ontology Biological Process (GO:BP) terms using clusterProfiler with Benjamini-Hochberg correction (FDR < 0.05). Twenty GO:BP terms were significantly enriched (Figure 5). The most strongly enriched term was transmembrane transport (GeneRatio ≈ 0.12, p.adjust < 1×10⁻⁹), which also had the largest gene count, indicating broad remodelling of membrane transport capacity in mature velum cells. Additional highly enriched terms included generation of precursor metabolites and energy, energy derivation by oxidation of organic compounds, cellular respiration, and aerobic respiration, all with adjusted p-values below 1×10⁻⁸. This pattern strongly supports a metabolic shift toward oxidative phosphorylation and respiratory energy generation in the mature biofilm, consistent with the aerobic lifestyle that characterises established flor velum (Mardanov et al., 2020).
At the lower end of enrichment significance, terms related to fermentative and glycolytic metabolism were also over-represented, including fermentation, hexose catabolic process, glucose catabolic process, and glycolytic process through fructose-6-phosphate and glucose-6-phosphate. The co-enrichment of both respiratory and fermentative pathways among the DEG set reflects the metabolic complexity of velum maturation, where genes in both pathways are simultaneously regulated, likely representing a transition from fermentation-dominant to respiration-dominant metabolism. Taken together, the ORA results indicate that velum maturation in S. cerevisiae flor yeast involves coordinated transcriptional induction of oxidative energy metabolism, ion homeostasis, and transmembrane transport, processes essential for survival and function at the air-liquid interface.

<img width="2400" height="2100" alt="GO_ORA_mature_vs_early" src="https://github.com/user-attachments/assets/878fc6f0-6917-492e-a229-313cee0ebf87" />
Figure 5. GO Biological Process enrichment dot plot (Mature vs Early, ORA). The x-axis shows the GeneRatio (proportion of DEGs annotated to each term); dot size is proportional to the number of DEGs in the term; dot colour indicates adjusted p-value (red = more significant, blue = less significant). Only terms with FDR < 0.05 are shown. Analysis performed using clusterProfiler with org.Sc.sgd.db and BH correction.

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
