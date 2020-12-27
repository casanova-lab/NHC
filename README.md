# NHC (Network-based Heterogeneity Clustering)
NHC: A computational approach to detect physiological homogeneity in the midst of genetic heterogeneity

## Background
The human genetic dissection of a growing range of clinical phenotypes is facing the challenge of genetic heterogeneity. Emerging data suggest that physiological homogeneity connects the products of genes whose variations underlie a given phenotype. While gene burden approaches have been used to identify genetic signals in case-control studies, it is underpowered in genetically heterogenous cohorts, and does not account for biological relevance anong genes. We therefore developed a computational genome-wide method to detect physiological homogeneity in the midst of genetic heterogeneity. Simulation studies and real patient cohort studies showed that our Network-based Heterogenity Clustering (NHC) method is able to systematically converge genes of biological proximity on a background biological interaction network, and capture gene clusters that harbor presumably deleterious mutations, in an efficient and unbiased manner. 

## Method Introduction
We are providing the codes for gene clustering for (1) patient cohort only and (2) patient cohort vs control cohort.

Although the goal of our method is to detect presumably deleterious mutations in genes that are clustered with close biological relevance from a cohort of patients, it is difficult to provide the code starting from varaint-level processing, as the variant data format and the variant filtration criteria vary hugely from one to another. Therefore, we leave the variant-level processing to each user, and the users need to prepare the gene list for all the individuals studied. Our code works on gene-level, and converges the gene candidates into gene clusters with pathway enrichment.

### Flowchart
![Image of NHC_Fig1](http://shiva.rockefeller.edu/NHC/NHC_GitHub_Fig_1.png)

### Brief Description
- A large-scale network of biological interaction is established, by integrating the human protein-protein interations (PPIs) from BioGRID, IntAct and REACTOME databases. PPIs are weighted by using the scores from STRING database to represent the level of biological relevance among genes. We obtained an edge-weighted background biological network of 202,057 PPIs for 15,585 human genes.

- By providing the gene list in patient cohort after variant filtration, our method traverses all genes of all patients in the edge-weighted background network, and iteratively converge genes with biological proximity into gene clusters. The algorithm starts from one gene of one patient, and iteratively searches for the closest gene in the rest of the patients that is above the edge-weight cutoff, where we used a stringent cutoff (STRING score ≥ 0.99 as default) to converge the gene clusters of the highest biological relevance. Each round of clustering stops when all patients have been visited or no other gene in the unvisited patients is above the edge-weight cutoff, and then outputs one gene cluster and its corresponding patient cluster. The algorithm will resume the clustering by starting from another gene of this patient, until every gene of every patient has been used as the start point for clustering.

- Converged gene clusters then iteratively merges two clusters, if one is a superset/subset of another, or two most-overlapping clusters sharing over 50% (default) genes, thus to output a reduced number of gene clusters that are more distinct to each other.

- Determine the statistical significance of each gene cluster in patients by principle components (PC) adjusted cluster-level enrichment versus controls. *(this step will be skipped, if users choose to run our method for patient cohort only)*

- Enrichement will be conducted on 1,720 pathways (187 KEGG pathways and 1,533 REACTOME pathways), collected from MSigDB database. We use pvalue 1e-5 as the significance cutoff for pathway enrichment, and assign the most-enriched pathway (the pathway with the lowest pvalue) to surrogate the primary physiological nature of each gene cluster.

## Usage
### Dependency
The code is written in python3, requiring python packages scipy *(both)* and *rpy2 (only in patient_vs_control)*

### Illustration  
![Image of NHC_Fig2](http://shiva.rockefeller.edu/NHC/NHC_GitHub_Fig_2.png)

### File Format
**Input:** Gene list in patients and controls *(example: test_patients.txt, test_controls.txt)*
- tab-delimited text file, including header line
- column 1: sample ID
- column 2: gene list separated by ',' without space

**Input:** PC (principal component) table for patients and controls *(example: test_pc.txt)*
- tab-delimited text file, including header line
- column 1: sample ID
- column 2-4: first 3 PC values for each sample

**Output:** Gene clusters in patients *(example: test_patients_output_gene_clusters.patient_only.txt, test_patients_output_gene_clusters.patient_vs_controls.txt)*
- tab-delimited text file, including header line
- columns:
  - Cluster ID
  - Number of Genes
  - Number of Patients
  - Gene Cluster
  - Patient Cluster
  - Cluster pValue *(only in patient_vs_control)*
  - Number of Pathways
  - Pathway List
  - Top Pathway
  - Top Pathway pValue

### Commands
**Default parameters:**
```
python NHC_code_patient_only.py -p test_patients.txt
```
```
python NHC_code_patient_vs_controls.py -p test_patients.txt -c test_controls.txt -pc test_pc.txt
```

**Customizable parameters:**
```
python NHC_code_patient_only.py -p <txt> -w <txt> -b <int> -m <float>
```
```
python NHC_code_patient_vs_controls.py -p <txt> -c <txt> -pc <txt> -w <float> -b <int> -m <float>
```

### Parameters
Parameter | Type | Description | Default
----------|------|-------------|--------------
*-p*|text file|gene list per ***p***atient (incl. header line)|na
*-c*|text file|gene list per ***c***ontrol (incl. header line)|na
*-pc*|text file|3 ***pc*** value for all samples (incl. header line)|na
*-w*|float|edge-***w***eight cutoff, based on STRING score [0~1]|0.99
*-b*|int|remove hu***b*** genes with high connectivity, use 0 to include all genes|50
*-m*|float|***m***erge overlapped (common genes/union genes) gene clusters|0.5

***Note:***
- *Strigent edge-weight cutoff (defalut 0.99) is to converge the gene clusters of the highest biological relevance. If the patient cohort is small or the gene candidates are few, the users could relax the edge-weight cutoff to 0.95, 0.9, but no lower than 0.7 (as STRING determines 0.7 as low-confidence cutoff).*
- *Hub gene removal is to avoid super-huge clusters that are converged due to the hub genes have huge amount of interacting genes. The connectivity of each gene is determined by the number of PPIs above STRING score 0.9 (NHC_data_connectivity.txt). The default value (-b 50) means: we are skipping the genes having more than 50 PPIs with edge-weight > 0.9 for clustering. If users want to include all genes for clustering, use (-b 0).*

## References
- *Zhang P. et al.* A computational approach to detect physiological homogeneity in the midst of genetic heterogeneity. (2020)
- *Casanova J.L. & Abel L.* The human genetic determinism of life-threatening infectious diseases: genetic heterogeneity and physiological homogeneity? *Hum Genet* (2020) [PubMed](https://pubmed.ncbi.nlm.nih.gov/32462426/)
- *McClellan J. & King M.C.* Genetic heterogeneity in human disease. *Cell* (2010) [PubMed](https://pubmed.ncbi.nlm.nih.gov/20403315/)
- *Povysil G. et al.* Rare-variant collapsing analyses for complex traits: guidelines and applications. *Nat Rev Genet* (2019) [PubMed](https://pubmed.ncbi.nlm.nih.gov/31605095/)
- *Itan, Y. et al.* The human gene connectome as a map of short cuts for morbid allele discovery. *PNAS* (2013) [PubMed](https://pubmed.ncbi.nlm.nih.gov/23509278/)

## Contact
> **Author:** Peng Zhang

> **Email:** pzhang@rockefeller.edu

> **Laboratory:** St. Giles Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA
