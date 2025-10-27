# SpLin v1.0
## 1. Outline
SpLin is an R package specifically designed for high-resolution spatially resolved transcriptomics, integrating three core functionalities to enable end-to-end analysis from tissue morphological correction to functional hotspot resolution. First, spatial data is partitioned into "windows"  using customized grid-based or fragment-based strategies to construct spatial array window cells, achieving precise mapping and co-localization of gene modules with spatial coordinates and tissue structures. Next, identification of HGMSRs based on the Integrated Genomic Score Index (IGSI) enables genome-wide detection of key functional domains exhibiting high expression variability and enrichment. Finally, dynamic hotspot tracking—Module Co-activity Spot (MCS) analysis —captures the most active micro-regions within modules in real-time via energy gradients and co-activity scores, providing a high-resolution framework for deciphering spatiotemporal regulatory mechanisms in complex biological processes. To investigate regional features of the intestin, we also developed an unrolling algorithm that converts “Swiss-roll” intestinal tissue into a linear layout, operating on both the expression matrix of the rolled sample and the corresponding ssDNA or H&E-stained images. 
## 2. Usage
`cd SpLin\documentation and see the SpLin-document`
## 3. Workflow
![workflow1](https://github.com/scholarLW/SpLin/blob/main/documentation/workflow_AI_a_01.png)
![workflow1](https://github.com/scholarLW/SpLin/blob/main/documentation/workflow_AI_b_01.png)
