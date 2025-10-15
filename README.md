# DISCOdiagnostics
Scripts for: From Concept to Challenge: Microfluidic Approach for Cell-Based Non-Invasive Testing

Non-invasive prenatal testing (NIPT) using cell free fetal DNA has transformed prenatal care by allowing risk assessment for a number of fetal genetic conditions without posing risk to the pregnancy. However, its limitations, including low fetal fraction and high DNA fragmentation, restrict its capacity for comprehensive genomic analysis. As an alternative, cell-based NIPT offers the potential to capture intact fetal cells for more detailed genetic interrogation. Cervical mucosa is a promising source for fetal DNA as it can be collected non-invasively as early as 5 weeks into the pregnancy and contains intact fetal cells. 

We performed single cell mRNA squencing and SNP nanopore sequencing of transcervical swabs and placental samples from pregant donors at multiple gestational ages to identify potential fetal cells. 

The provided scripts go through mapping, parsing fasta headers (barcodes and UMIs), generating count data, making figures and performing differential gene expression (DGE):

## General Workflow
The provided scripts go through the initial processing and merging of the outputs from CellRanger in R using Seurat. Analysis of FastQ files were performed by the Princess Margaret Genomics Centre.

Initial loading, normalizing, annotation, and visualizing of individual datasets into a Seurat Object
Merging swab datasets and analyzing cellular profiles
Mergign placental and swab datasets and generation of fetal cell marker panel

## Publication
*In progress*
Dylan Siriwardena, Michael Dryden, M. Dean Chamberlain, Chloe Taylor, Louise Dupoiron, Farhana Abbas, Julian Lamanna, David Chitayat, Aaron R. Wheeler and Elena Greenfeld. From Concept to Challenge: Microfluidic Approach for Cell-Based Non-Invasive Testing. 2025.

## Raw data
Raw data is stored at Arrayexpress: 
