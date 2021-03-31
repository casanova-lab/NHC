# NHC (Network-based Heterogeneity Clustering)
NHC: A computational approach to detect physiological homogeneity in the midst of genetic heterogeneity *(LICENSE: CC BY-NC-ND 4.0)*

## Background
The human genetic dissection of a growing range of clinical phenotypes is facing the challenge of genetic heterogeneity. Emerging data suggest that physiological homogeneity connects the products of genes whose variations underlie a given phenotype. While gene burden approaches have been used to identify genetic signals in case-control studies, it is underpowered in genetically heterogenous cohorts, and does not account for biological relevance between genes.

We therefore developed a computational genome-wide method (NHC) to systematically converge genes of biological proximity on a background biological interaction network, and capture the gene clusters that harbor presumably deleterious mutations, in an efficient and unbiased manner. NHC method is suitable for rare and common diseases, that have an homogeneous clinical phenotype, and are likely caused by rare/uncommon variants with strong individual effects located in physiologically related genes, at least in a group of individuals with a given condition.

## Method Introduction
We are providing the codes for gene clustering for (1) case cohort only and (2) case-control cohort.

Although the goal of our method is to detect presumably deleterious mutations in genes with close biological relevance from a cohort of cases with the same disease, it is difficult to provide the code starting from variant-level processing, as the variant data format and the variant filtration criteria vary hugely from one lab to another, and from one study to another. Therefore, we leave the variant-level processing to the users, who need to prepare the gene list for all the individuals under study. Our code works on gene-level, and converges the genes carrying qualifying variants into gene clusters with pathway and gene ontology enrichment.

### Flowchart
![Image of NHC_Fig1](http://shiva.rockefeller.edu/NHC/NHC_GitHub_Fig1.png)

### Brief Description
- A large-scale network of human protein-protein interactions (PPIs) is established, based on BioGRID, IntAct and REACTOME databases. PPIs are weighted by using the scores from STRING database to represent the level of biological relevance between genes. We obtained an edge-weighted background biological network of 202,057 PPIs for 15,585 human genes.

- By providing the gene list in case cohort after variant filtration, NHC traverses all genes of all cases in the edge-weighted background network, and iteratively converge genes with biological proximity into gene clusters. The algorithm starts from one gene of one case, and iteratively searches for the closest gene in the rest of the cases that is above the edge-weight cutoff. We used STRING score ≥ 0.99, a stringent default, to converge the clusters of the highest biological relevance. Each round of clustering stops when all cases have been visited or no other gene in the unvisited cases is above the edge-weight cutoff, and then outputs one gene cluster and its corresponding case cluster. The algorithm will resume the clustering by starting from another gene of this case, until every gene of every case has been used as the start point for clustering once.

- The initial outputted gene clusters are then iteratively merged, if one is a superset/subset of another, or the two most-overlapping clusters sharing over 50% (default) genes, thereby reducing the number of gene clusters that are more distinct from each other.

- Determine the statistical significance of each gene cluster in cases versus controls by principle components (PC) adjusted cluster-level enrichment, which aims to reduce the effect of ethnic diversity. *(this step will be skipped, if users choose to run our method for case cohort only)*

- Enrichment will be conducted on pathways (187 KEGG and 1,533 REACTOME pathways) and gene ontologies (GO) (8,992 biological process (BP) and 2,812 molecular function (MF)). We use p-value 1e-5 as the significance cutoff for pathway/GO enrichment.

- In order to deal with large cohorts of cases, we also provide a boost version (NHC-boost) of the gene clustering algorithm, which follows the same concept of the original algorithm, but traverses each gene of a specific case only once. In other words, if a given gene of a specific case has been clustered into one cluster, it will not be traversed and clustered again in the rest of clustering iterations. The performance of NHC-boost may mildly decrease, but significantly increases the computation efficiency.

## Usage
### Dependency
The code is written in python3, requiring python packages scipy *(case only, case-control)* and *rpy2 (case-control)*

### Illustration  
![Image of NHC_Fig2](http://shiva.rockefeller.edu/NHC/NHC_GitHub_Fig2.png)

### File Format
**Input:** Gene list in cases and controls *(example: test_cases.txt, test_controls.txt)*
- tab-delimited text file, including a header line
- column 1: sample ID
- column 2: gene list separated by ',' without space

**Input:** PC (principal component) table for cases and controls *(example: test_pc.txt)*
- tab-delimited text file, including a header line
- column 1: sample ID
- column 2-4: first 3 PC values for each sample *(if PCs are unavailable, use 1 for all, assuming no ethnic diversity)*

**Output:** Gene clusters converged in cases *(example: output_case_only.txt, output_case_control.txt)*
- tab-delimited text file, including a header line
- columns:
  - Cluster ID
  - Number of Genes
  - Number of Cases
  - Gene Cluster
  - Case Cluster
  - Cluster pvalue *(only in patient_control)*
  - Number of Pathways
  - Pathway List
  - Top Pathway
  - Top Pathway pvalue
  - GO BP List
  - GO MF List

### Commands
**Default parameters:**
```
python NHC_case_only.py -case test_cases.txt
```
```
python NHC_case_control.py -case test_cases.txt -ctl test_contorls.txt -pc test_pc.txt
```

**Customizable parameters:**
```
python NHC_case_only.py -case <txt> -w <float> -b <int> -m <float> -o <txt>
```
```
python NHC_case_control.py -case <txt> -ctl <txt> -pc <txt> -w <float> -b <int> -m <float> -o <txt>
```

### Parameters
Parameter | Type | Description | Default
----------|------|-------------|--------------
*-case*|file|gene list per case (incl. a header line)|na
*-ctl*|file|gene list per control (incl. a header line)|na
*-pc*|file|3 pc value for all samples (incl. a header line)|na
*-w*|float|edge-weight cutoff, based on STRING score [0~1]|0.99
*-b*|int|remove hub genes with high connectivity, use 0 to include all genes|50
*-m*|float|merge overlapped clusters (overlapping ratio = common/union genes)|0.5
*-o*|text|output filename|output_timestamp.txt

***Note:***
- *Stringent edge-weight cutoff (default 0.99) is used to converge the gene clusters of the highest biological relevance. If the case cohort is small or the gene candidates are few, the users could relax the edge-weight cutoff to 0.95, 0.9, but no lower than 0.7 (as STRING determines 0.7 as low-confidence cutoff).*
- *Hub gene removal is to avoid giant clusters that are formed due to the hub genes have large amount of interacting genes. The connectivity of each gene is determined by the number of PPIs above STRING score 0.9 (Data_connectivity.txt). The default value (-b 50) means: skipping the genes having more than 50 PPIs with edge-weight>0.9 for clustering. If users want to include all genes for clustering, use (-b 0).*
- *NHC-boost has the same setting for parameters and the same output format, just call NHC-boost_case_only.py or NHC-boost_case_control.py.*

## References
- *Zhang P. et al.* A computational approach to detect physiological homogeneity in the midst of genetic heterogeneity. (2021)
- *Casanova J.L. & Abel L.* The human genetic determinism of life-threatening infectious diseases: genetic heterogeneity and physiological homogeneity? *Hum Genet* (2020) [PubMed](https://pubmed.ncbi.nlm.nih.gov/32462426/)
- *McClellan J. & King M.C.* Genetic heterogeneity in human disease. *Cell* (2010) [PubMed](https://pubmed.ncbi.nlm.nih.gov/20403315/)
- *Povysil G. et al.* Rare-variant collapsing analyses for complex traits: guidelines and applications. *Nat Rev Genet* (2019) [PubMed](https://pubmed.ncbi.nlm.nih.gov/31605095/)
- *Itan Y. et al.* The human gene connectome as a map of short cuts for morbid allele discovery. *PNAS* (2013) [PubMed](https://pubmed.ncbi.nlm.nih.gov/23509278/)

## Contact
> **Author:** Peng Zhang, Ph.D.

> **Email:** pzhang@rockefeller.edu

> **Laboratory:** St. Giles Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA
