---
layout: default
title: Project Superdarks
---
<!-- stager: preserve -->

# Project Superdarks

This project finds structural homologs ("darks") for a query protein across the entire AlphaFold Database using distributed computing, then enriches and networks the hits for discovery.

- [code](code/)
- [data](data/)

## 🌌 What Superdarks Does

**Superdarks** is a large-scale **structural homology discovery platform**. Starting from a single query PDB structure (for example, a 7-transmembrane receptor), it performs pairwise **TM-align** of the query against **214,528,851 AlphaFold Database predictions**, distributes the work across up to **1,000 HPC compute nodes**, and post-processes the hits through a four-stage informatics pipeline that produces **ranked hit lists, UpSet diagrams, contour plots, UniProt + InterPro annotations, BLAST-based sequence networks, and Cytoscape-ready subgraphs.**

## 🔍 Key Capabilities

- **Distributed Structural Alignment**  
  Submits per-node **TM-align** jobs to the Triton or Pegasus clusters (or BU SCC) via LSF `bsub`, with a job parser that consolidates hits and supports resubmission of failed jobs.

- **Filtering by Coverage and TM-Score**  
  Parses `rank--coverage-tm--UNIPROT-…-full.pdb` filenames and applies configurable `min_tm_score` and `min_residues` thresholds to distinguish meaningful structural matches from background noise.

- **Hit Consolidation and Re-alignment**  
  Merges per-node result objects, then optionally re-aligns top hits to produce aligned PDBs suitable for figure generation.

- **UpSet Diagram Analysis**  
  Builds UpSet intersection diagrams of domain membership for both **filtered** and **unfiltered** hits, producing publication-ready PNG / SVG output and CSV intersection tables.

- **Score Contour Plotting**  
  Headless plotting pipeline that emits **2D and 3D scatter plots** of coverage vs TM-score vs rank, density contours, and a companion Excel table — optionally per-rank.

- **UniProt + InterPro Enrichment**  
  Streams Swiss-Prot and TrEMBL `.dat.gz` entries (via `pigz` / `gzip`), merges them with hit listings and InterPro annotations, and emits a unified TSV with gene, organism, superkingdom, and protein-function columns.

- **7TM / GPCR Filtering**  
  Pandas-based TSV filter with flexible 7TM+GPCR keyword matching (including a GCR-family rule) that writes per-match Excel workbooks up to the 1,048,576-row sheet limit.

- **Sequence Chunking + BLASTP Networks**  
  Filters FASTAs by UniProt-ID list, chunks sequences into pair-FASTAs, submits BLASTP jobs via LSF, then builds a **gzipped nodes/edges network** (FAST in-process or TURBO/cluster) ready for Cytoscape.

- **Go-Forward Network Traversal**  
  Applies the Isom-Lab **goFo** recursive walk to the thresholded BLAST similarity graph, with filters by percent identity, e-value, or bitscore.

- **Conservation-to-Structure Mapping**  
  Performs Clustal-Omega pairwise traceback to map sequence conservation onto a reference PDB, enabling 3D visualization of evolutionary signal.

## 🛠️ Additional Functionalities

- **Foldseek Wrapper (GUI)**  
  `foldseek_individual_query_gui.py` is a Tkinter front-end for running `foldseek createdb` and search, with stale-index detection and macOS-safe file filters.

- **fpocket Wrapper**  
  `fpocket_python_wrapper.py` is a thin Python wrapper around the `fpocket` cavity-detection binary that gathers `*_out/` directories into a clean result set.

- **Cluster-Agnostic Configuration**  
  Each script's CONFIG block (or `getopt` flags) lets you swap between Triton, Pegasus, BU SCC, or a local workstation by changing one or two variables at the top of the file.

## 📚 Citations

Superdarks is a pipeline built on several foundational methods. If you use it in your work, please cite the underlying tools along with any Isom-Lab publications describing the overall workflow (in preparation).

1. Zhang Y, Skolnick J. *TM-align: a protein structure alignment algorithm based on the TM-score.* Nucleic Acids Res. 2005 Apr 22;33(7):2302-9. doi: [10.1093/nar/gki524](https://doi.org/10.1093/nar/gki524). PMID: [15849316](https://pubmed.ncbi.nlm.nih.gov/15849316)
2. Jumper J, Evans R, Pritzel A, et al. *Highly accurate protein structure prediction with AlphaFold.* Nature. 2021 Aug;596(7873):583-589. doi: [10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2). PMID: [34265844](https://pubmed.ncbi.nlm.nih.gov/34265844)
3. Varadi M, Anyango S, Deshpande M, et al. *AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models.* Nucleic Acids Res. 2022 Jan 7;50(D1):D439-D444. doi: [10.1093/nar/gkab1061](https://doi.org/10.1093/nar/gkab1061). PMID: [34791371](https://pubmed.ncbi.nlm.nih.gov/34791371)
4. van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. *Fast and accurate protein structure search with Foldseek.* Nat Biotechnol. 2024 Feb;42(2):243-246. doi: [10.1038/s41587-023-01773-0](https://doi.org/10.1038/s41587-023-01773-0). PMID: [37156916](https://pubmed.ncbi.nlm.nih.gov/37156916)
5. Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ. *Gapped BLAST and PSI-BLAST: a new generation of protein database search programs.* Nucleic Acids Res. 1997 Sep 1;25(17):3389-402. doi: [10.1093/nar/25.17.3389](https://doi.org/10.1093/nar/25.17.3389). PMID: [9254694](https://pubmed.ncbi.nlm.nih.gov/9254694)
6. Sievers F, Higgins DG. *Clustal Omega for making accurate alignments of many protein sequences.* Protein Sci. 2018 Jan;27(1):135-145. doi: [10.1002/pro.3290](https://doi.org/10.1002/pro.3290). PMID: [28884485](https://pubmed.ncbi.nlm.nih.gov/28884485)
7. Le Guilloux V, Schmidtke P, Tuffery P. *Fpocket: an open source platform for ligand pocket detection.* BMC Bioinformatics. 2009 Jun 2;10:168. doi: [10.1186/1471-2105-10-168](https://doi.org/10.1186/1471-2105-10-168). PMID: [19486540](https://pubmed.ncbi.nlm.nih.gov/19486540)
