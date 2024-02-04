# NHC (Network-based Heterogeneity Clustering)
NHC: A computational approach to detect physiological homogeneity in the midst of genetic heterogeneity *(LICENSE: CC BY-NC-ND 4.0)*

## Introduction
The human genetic dissection of a growing range of clinical phenotypes is facing the challenge of genetic heterogeneity. Emerging data suggest that physiological homogeneity connects the gene products whose variations underlie a given phenotype. Gene burden tests are to identify genetic signals in case-control studies by assuming genetic homogeneity, which can be underpowered in genetically heterogeneous cohorts.

We developed NHC method to systematically converge genes of biological proximity on a background protein-protein interaction network, and to capture the gene clusters that harbor presumably deleterious variants, in an unbiased manner. NHC method is suitable for studying the patient cohort with a homogeneous clinical phenotype, which is likely caused by rare or uncommon variants with strong individual effects in physiologically related genes.

### Flowchart
<img src="https://hgidsoft.rockefeller.edu/NHC/Figure_NHC_v3.png" width="77%" height="77%">

### Description
- A large-scale network of human protein-protein interactions (PPIs) is established, based on STRING, BioGRID and REACTOME databases. PPIs are required to be physical, and weighted by STRING scores to represent the quantitative biological relevance between genes. The edge-weighted background biological network includes 158,326 PPIs for 13,612 human genes.

- By providing the sample-gene-variant data from a case cohort after variant filtration, NHC traverses all genes of all cases in the edge-weighted background network, and iteratively converges genes with biological proximity into gene clusters. The algorithm starts from one gene of one case, and iteratively searches for the closest gene in the rest of the cases that is above the edge-weight cutoff (default: -edge 0.99). A stringent default is to converge the clusters of the highest biological relevance. Each round of clustering stops when all cases have been visited or no gene in the unvisited cases is above the edge-weight cutoff, and then outputs one gene cluster and its corresponding case cluster. The algorithm will resume the clustering by starting from the next gene of this case, until every gene of every case has been used as the starting point for clustering.

- The initial output gene clusters are then iteratively merged, if one is a superset/subset of another, or the two most-overlapping clusters share over 50% (default: -merge 0.5) genes, thereby generating gene clusters that are more distinct from each other.

- Determine the statistical significance of each gene cluster in a case-vs-control study, by principal components (PC) adjusted cluster-level enrichment.

- Geneset enrichment based on MSigDB Hallmark (50 genesets), KEGG Pathway (805 genesets), Reactome Pathway (1670 genesets), Wiki Pathway (791 genesets), GO Biological Process (7375 genesets), GO Molecular Function (1753 genesets). NHC uses adjusted p-value 1e-5 as the significance cutoff.

- To deal with the large number of samples and genes, NHC provides a boost version, which follows the same concept of the original clustering algorithm, but traverses each gene of a specific case only once. Its clustering is not as greedy as the normal version, and its performance may mildly decrease, but it significantly increases the computation efficiency.

## News
- 02/2024: NHC official version-3 was released, with new features: accepting variant-level input, outputting network files for visualization, integrated case-only and case-vs-control modes; integrated normal and boost versions; supported more geneset enrichment; and updated background protein-protein interaction network.
- 04/2023: NHC official version-2 was released, with new features: updated background protein-protein interaction network; and updated cluster-level enrichment test.
- 06/2021: "*A computational approach for detecting physiological homogeneity in the midst of genetic heterogeneity*" that introduces NHC method was published in [*The American Journal of Human Genetics*](https://www.cell.com/ajhg/fulltext/S0002-9297(21)00154-3).
- 12/2020: NHC official version-1 was released.
- 07/2020: NHC prototype was developed.

## Usage
Current version: version-3
### Dependency
The code is written in python3, requiring python packages [*scipy*](https://scipy.org/install/) and [*rpy2*](https://rpy2.github.io/doc/latest/html/index.html).

### Input / Output
**Input:** 

- The sample-gene-variant data under study (example: *test_input.txt*)
  - Tab-delimited text file, including a header line
  - Column 1: group (case or control)
  - Column 2: sample ID
  - Column 3: gene
  - Column 4 onwards: can be any variant-level information, which will be entirely extracted in the output. 
- The principal component (PC) table for all samples (example: *test_pc.txt*)
  - Tab-delimited text file, including a header line
  - Column 1: sample ID
  - Column 2-4: first three PCs for each sample (if no PC, use 1 for all, assuming no ethnic diversity)

**Output:**

Each run will create a new folder in the given path, with the folder name (*NHC_output_timestamp_suffix*)
- *NHC_input_parameters.txt*
  - A record of the parameters used in this run
- *NHC_output_gene_clusters.txt*
  - The final gene clusters converged from the input data, with the following columns:
  - Cluster ID
  - Gene Count, Gene Cluster
  - Case Count, Case Cluster
  - Cluster pValue
  - Geneset Enrichment (MSigDB_Hallmark, KEGG_Pathway, Reactome_Pathway, Wiki_Pathway, GO_BiologicalProcess, GO_MolecularFunction)
- Folder *variant_files*
  - *NHC_output_gene_cluster_1...n_variants.txt* (provides the sample-gene-variant data for each gene cluster, following the exact format of the input file)
- Folder *network_files*
  - *NHC_output_gene_cluster_1...n_network.csv* (provides the network file for each gene cluster, that can be imported to Gephi or Cytoscape software for visualization)
  - *NHC_output_gene_cluster_1...n_node.csv*, (provides the case count and variant count per gene in each gene cluster, that can be used to visualize the gene node size)

### Command
**Command Parameters:**
```
python NHC.py -path -input -pc -mode -edge -hub -merge -boost -network -suffix
```
**Example Command:**
```
python NHC.py -path /x/y/z/ -input test_intput.txt -pc test_pc.txt -mode 2 -edge 0.99 -hub 100 -merge 0.5 -boost N -network Y -suffix test
```
Parameter | Type | Description | Default
----------|------|-------------|--------------
*-path*|text|absolute path of the input data|na
*-input*|file|input file for samples, genes, and variants (including header)|na
*-pc*|file|three principal components for all samples (including header)|na
*-mode*|int|1 for case-only analysis, 2 for case-vs-control analysis|1
*-edge*|float|edge weight cutoff, range: 0.7~1|0.99
*-hub*|int|remove hub genes with high connectivity, use 0 to keep all genes|100
*-merge*|float|merge overlapped gene clusters, range: 0~1|0.5
*-boost*|text|Y or N to use boost version|N
*-network*|text|Y or N to generate network files for visualization|N
*-suffix*|text|suffix for output folder|na

***Note:***
- *Stringent edge-weight cutoff (default: -edge 0.99) is used to converge the gene clusters of the highest biological relevance. If the case cohort is small or the gene candidates are few, then users could relax the edge-weight cutoff to 0.95 or 0.9, but no lower than 0.7.*
- *Hub gene removal is to avoid generating giant clusters due to the large number of interactions with hub genes. The connectivity of each gene (default: -hub 100) is determined by the number of its high-confidence PPIs, meaning to skip the genes having more than 100 PPIs for clustering. If users want to include all genes for clustering, use (-hub 0).*
- *Boost version has the same input/output format, with faster computing speed, especially when the input data (number of samples and genes) is huge. Its clustering algorithm is not as greedy as the normal version.*
- *Network option provides the network files and gene node size files, that can be imported to Gephi or Cytoscape software, for visualizing each cluster.*

## References
- *Zhang P. et al.* A computational approach to detect physiological homogeneity in the midst of genetic heterogeneity. [*Am J Hum Genet* (2021)](https://www.cell.com/ajhg/fulltext/S0002-9297(21)00154-3)
- *Casanova J.L. & Abel L.* The human genetic determinism of life-threatening infectious diseases: genetic heterogeneity and physiological homogeneity? [*Hum Genet* (2020)](https://pubmed.ncbi.nlm.nih.gov/32462426/)
- *McClellan J. & King M.C.* Genetic heterogeneity in human disease. [*Cell* (2010)](https://pubmed.ncbi.nlm.nih.gov/20403315/)
- *Povysil G. et al.* Rare-variant collapsing analyses for complex traits: guidelines and applications. [*Nat Rev Genet* (2019)](https://pubmed.ncbi.nlm.nih.gov/31605095/)
- *Itan Y. et al.* The human gene connectome as a map of shortcuts for morbid allele discovery. [*PNAS* (2013)](https://pubmed.ncbi.nlm.nih.gov/23509278/)

## Contact
> **Author:** Peng Zhang, Ph.D.

> **Email:** pzhang@rockefeller.edu

> **Laboratory:** St. Giles Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA
