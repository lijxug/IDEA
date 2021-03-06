[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](https://github.com/lijxug/IDEA/blob/master/LICENSE)

# IDEA

Single cell RNA-seq has enabled high-resolution characterization of molecular signatures of tumor-infiltrating lymphocytes. However, analyses at the transcript isoform level are rarely reported. As alternative splicing is critical to T cell differentiation and activation, here we proposed a computational method named as IDEA to comprehensively detect and annotate differentially used isoforms across cell subtypes.
This repository works as an example to reproduce figures in the [article](https://doi.org/10.1101/2020.01.29.924308). The demo data, which contains the expression data of CD8+ T cells and transcripts of genes presented in the paper, is available in `/demo/`.

Run the codes in the following order:

- [01.calculateEnrichment.Rmd](01.calculateEnrichment.Rmd) perform the enrichment analysis 
- [02.Annotation_and_Event_detection.Rmd](02.Annotation_and_Event_detection.Rmd) add functional and structural annotations to transcripts
- [03.plotEnrichmentofGenes.Rmd](03.plotEnrichmentofGenes.Rmd) drawing the organized plot. A successful run results in plots like the one shown below.

![CD8T_CD8A](CD8T_CD8A.png)


