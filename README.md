# NHC (Network-based Heterogeneity Clustering)
NHC: A computational approach to detect physiological homogeneity in the midst of genetic heterogeneity *(LICENSE: CC BY-NC-ND 4.0)*

## Background
The human genetic dissection of a growing range of clinical phenotypes is facing the challenge of genetic heterogeneity. Emerging data suggest that physiological homogeneity connects the gene products whose variations underlie a given phenotype. Gene burden tests are to identify genetic signals in case-control studies by assuming genetic homogeneity, which are underpowered in genetically heterogenous cohorts.

We developed NHC method to systematically converge genes of biological proximity on a background biological interaction network, and to capture the gene clusters that harbor presumably deleterious variants, in an unbiased manner. NHC method is suitable for studying the patient cohort with a homogeneous clinical phenotype, which are likely caused by rare/uncommon variants with strong individual effects in physiologically related genes.

## Introduction
We are providing the codes for gene clustering for (1) case cohort vs control cohort, and (2) case cohort only.

Although the goal of our method is to detect presumably deleterious mutations in genes with close biological relevance from a case cohort with the same disease, it is difficult to provide the code starting from variant-level processing, as the variant data format and the variant filtration criteria vary hugely from one lab to another, and from one study to another. Therefore, we leave the variant-level processing to the users, who need to prepare the candidate gene list for all the individuals under study. Our code works on gene-level, and converges the genes carrying the qualified variants into gene clusters with pathway and gene ontology enrichment.

### Flowchart
<img src="https://hgidsoft.rockefeller.edu/NHC/Figure_NHC_1.png" width="60%" height="60%">

### Description
- A large-scale network of human protein-protein interactions (PPIs) is established, based on STRING, BioGRID and REACTOME databases. PPIs are required to be physical, and then weighted by STRING scores to represent the biological relevance between genes. We built an edge-weighted background biological network of 157,205 PPIs for 13,283 human genes.

- By providing the candidate gene list in case cohort after variant filtration, NHC traverses all genes of all cases in the edge-weighted background network, and iteratively converge genes with biological proximity into gene clusters. The algorithm starts from one gene of one case, and iteratively searches for the closest gene in the rest of cases that is above the edge-weight cutoff (default: 0.99). A stringent default is to converge the clusters of the highest biological relevance. Each round of clustering stops when all cases have been visited or no gene in the unvisited cases is above the edge-weight cutoff, and then outputs one gene cluster and its corresponding case cluster. The algorithm will resume the clustering by starting from the next gene of this case, until every gene of every case has been used as the start point for clustering.

- The initial output gene clusters are then iteratively merged, if one is a superset/subset of another, or the two most-overlapping clusters sharing over 50% (default) genes, thereby generating gene clusters that are more distinct from each other.

- Determine the statistical significance of each gene cluster in cases versus controls by principal components (PC) adjusted cluster-level enrichment. *(this step will be skipped, if users run case_only)*

- Enrichment on pathways (187 KEGG and 1,533 REACTOME pathways) and gene ontologies (GO) (8,992 biological process (BP) and 2,812 molecular function (MF)). We use p-value 1e-5 as the significance cutoff.

- In order to deal with large number of cases, we also provide a boost version (NHCboost), which follows the same concept of the original algorithm, but traverses each gene of a specific case only once. In other words, if a given gene of a specific case has been clustered into one cluster, it will not be traversed and clustered again in the rest of clustering iterations. The performance of NHCboost may mildly decrease, but significantly increases the computation efficiency.

## News
- 04/2023: NHC official version-2 was released, with the updated background network and the updated enrichment test for gene clusters.
- 06/2021: "A computational approach for detecting physiological homogeneity in the midst of genetic heterogeneity" that introduces NHC method was published in [*AJHG*](https://www.cell.com/ajhg/fulltext/S0002-9297(21)00154-3).
- 12/2020: NHC official version-1 was released.

## Usage
Current version: version-2
### Dependency
The code is written in python3, requiring python packages scipy *(case_only, case_control)* and *rpy2 (case_control)*

### Illustration  
<img src="https://hgidsoft.rockefeller.edu/NHC/Figure_NHC_2.png" width="60%" height="60%">

### File Format
**Input:** Candidate gene list in cases and controls *(example: test_cases.txt, test_controls.txt)*
- tab-delimited text file, including a header line
- column 1: sample ID
- column 2: gene list separated by ',' without space

**Input:** PC (principal component) table for cases and controls *(example: test_pc.txt)*
- tab-delimited text file, including a header line
- column 1: sample ID
- column 2-4: first 3 PCs for each sample *(if no PC, use 1 for all, assuming no ethnic diversity)*

**Output:** Gene clusters converged in cases *(example: output_case_only.txt, output_case_control.txt)*
- tab-delimited text file, including a header line
- columns:
  - cluster ID
  - number of genes
  - number of cases
  - gene cluster
  - case cluster
  - cluster pvalue *(only in case_control)*
  - number of pathways
  - pathway list
  - top pathway
  - top pathway pvalue
  - GO BP list
  - GO MF list

### Commands & Parameters
**Default parameters:**
```
python NHC_case_control.py -case test_cases.txt -ctl test_controls.txt -pc test_pc.txt
```
```
python NHC_case_only.py -case test_cases.txt
```

**Customizable parameters:**
```
python NHC_case_control.py -case -ctl -pc -w -b -m -o
```
```
python NHC_case_only.py -case -w -b -m -o
```
Parameter | Type | Description | Default
----------|------|-------------|--------------
*-case*|file|gene list per case (incl. a header line)|na
*-ctl*|file|gene list per control (incl. a header line)|na
*-pc*|file|3 pc value for all samples (incl. a header line)|na
*-w*|float|edge-weight cutoff, based on STRING score [0.7~1]|0.99
*-b*|int|remove hub genes with high connectivity, use 0 to include all genes|100
*-m*|float|merge overlapped clusters (overlapping ratio = common/union genes)|0.5
*-o*|text|output filename|output_timestamp

***Note:***
- *Stringent edge-weight cutoff (default: 0.99) is used to converge the gene clusters of the highest biological relevance. If the case cohort is small or the gene candidates are few, then users could relax the edge-weight cutoff to 0.95 or 0.9, but no lower than 0.7 (as STRING determines 0.7 as confidence cutoff).*
- *Hub gene removal is to avoid giant clusters that are formed due to the large number of interactions with hub genes. The connectivity of each gene is determined by the number of PPIs above STRING score 0.9 (Data_Network_Connectivity.txt). The default value (-h 100) means: skipping the genes having more than 100 PPIs with edge-weight>0.9 for clustering. If users want to include all genes for clustering, use (-h 0).*
- *NHCboost has the same input/output format and the same parameter configurations, just call NHCboost_case_only.py or NHCboost_case_control.py.*

## References
- *Zhang P. et al.* A computational approach to detect physiological homogeneity in the midst of genetic heterogeneity. [*Am J Hum Genet* (2021)](https://www.cell.com/ajhg/fulltext/S0002-9297(21)00154-3)
- *Casanova J.L. & Abel L.* The human genetic determinism of life-threatening infectious diseases: genetic heterogeneity and physiological homogeneity? [*Hum Genet* (2020)](https://pubmed.ncbi.nlm.nih.gov/32462426/)
- *McClellan J. & King M.C.* Genetic heterogeneity in human disease. [*Cell* (2010)](https://pubmed.ncbi.nlm.nih.gov/20403315/)
- *Povysil G. et al.* Rare-variant collapsing analyses for complex traits: guidelines and applications. [*Nat Rev Genet* (2019)](https://pubmed.ncbi.nlm.nih.gov/31605095/)
- *Itan Y. et al.* The human gene connectome as a map of short cuts for morbid allele discovery. [*PNAS* (2013)](https://pubmed.ncbi.nlm.nih.gov/23509278/)

## Contact
> **Author:** Peng Zhang, Ph.D.

> **Email:** pzhang@rockefeller.edu

> **Laboratory:** St. Giles Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA
