# **Transcriptomic Profiling of *Saccharomyces cerevisiae* During Velum Biofilm Formation**
Author: Vian Lelo
Date created: January 20th, 2026, Last Updated: February 8th, 2026

# **1.0 Introduction**
## **1.1 Biological Background**
## **1.2 Objectives**
## **1.3 Methods Overview & Justification**

# **2.0 Methods**
## **2.1 Environment Setup**
```
```
### **2.1.1 Conda Environment**
```
```
### **2.1.2 Activation**
```
```
## **2.2 Data & Metadata**
```
```
### **2.2.1 Dataset**
```
```
### **2.2.2 Metadata File**
```
```
### **2.2.3 Data Download**
Reference Files:
- cdna.all.fa.gz: Transcriptome FASTA. Required for Salmon quantification; contains spliced transcript sequences where reads are actually mapped.
- dna.toplevel.fa.gz: Whole genome FASTA. Optional but recommended; used to build a decoy-aware Salmon index to reduce false mappings.
- .gtf.gz: Gene annotation (GTF). Required for tximport/DESeq2 to map transcript IDs → gene IDs for differential expression and GO/KEGG enrichment.

Created sample sheet for looping:
```

```
## **2.3 Quality Control**
```
```
### **2.3.1 FastQC & MultiQC**
```
```
### **2.3.2 Optional Trimming**
```
```
## **2.4 Reference & Indexing**
```
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
