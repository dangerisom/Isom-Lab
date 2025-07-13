---
layout: default
title: Project pHinder
---

# Project pHinder

This project explores structure-based relationsips using computational geometry.

- [code](code/)
- [data](data/)

## 🧬 What pHinder Does

**pHinder** is a computational toolkit for analyzing the 3D structure of proteins, with a focus on **side chain topology and burial**. Originally developed to identify **buried ionizable networks**, pHinder has evolved into a **general-purpose platform** for side chain classification and surface-based structural analysis.

## 🔍 Key Capabilities

- **Topological Classification of Side Chains**  
  Assigns each side chain as **core**, **margin**, or **exposed** based on its spatial relationship to a calculated molecular surface.

- **Molecular Surface Construction**  
  Builds and refines a triangulated surface around the protein using Cα-based Delaunay triangulation, enabling depth-based residue analysis.

- **Network Detection and Scoring**  
  Identifies contiguous **networks of side chains** (especially ionizable ones) based on proximity and burial, useful for exploring electrostatic or functional microenvironments.

- **Flexible Residue Sets**  
  Supports predefined or custom residue sets, including **ionizable**, **acidic**, **basic**, **polar**, **apolar**, or **user-defined** residues.

- **Burial-Based Filtering for Functional Inference**  
  Useful for inferring sites of potential structural stabilization, altered pKa behavior, or ligand interaction hotspots.

- **Generalization Beyond Electrostatics**  
  While originally tailored to **ionizable side chains**, pHinder now supports **broad structural analysis of any side chain class**, making it a versatile tool in protein structural biology.

## 🛠️ Additional Functionalities

- **Protein-Protein Interface Detection**  
  Automatically identifies interaction surfaces between two or more protein chains, facilitating the study of multimeric assemblies and complex formation.

- **Ligand and Drug Binding Site Identification**  
  Predicts buried and surface-accessible cavities that may serve as **potential ligand or drug binding sites**, supporting virtual screening workflows.

## 📚 Citations

If you use **pHinder** in your work, please cite the following publications that describe the method and its applications:

1. Isom DG, Sridharan V, Baker R, Clement ST, Smalley DM, Dohlman HG. *Protons as second messenger regulators of G protein signaling.* Mol Cell. 2013 Aug 22;51(4):531-8. doi: [10.1016/j.molcel.2013.07.012](https://doi.org/10.1016/j.molcel.2013.07.012). PMID: [23954348](https://pubmed.ncbi.nlm.nih.gov/23954348)
2. Isom DG, Dohlman HG. *Buried ionizable networks are an ancient hallmark of G protein-coupled receptor activation.* Proc Natl Acad Sci U S A. 2015 May 5;112(18):5702-7. doi: [10.1073/pnas.1417888112](https://doi.org/10.1073/pnas.1417888112). PMID: [25902551](https://pubmed.ncbi.nlm.nih.gov/25902551)
3. Isom DG, Sridharan V, Dohlman HG. *Regulation of Ras Paralog Thermostability by Networks of Buried Ionizable Groups.* Biochemistry. 2016 Jan 26;55(3):534-42. doi: [10.1021/acs.biochem.5b00901](https://doi.org/10.1021/acs.biochem.5b00901). PMID: [26701741](https://pubmed.ncbi.nlm.nih.gov/26701741)
4. Isom DG, Page SC, Collins LB, Kapolka NJ, Taghon GJ, Dohlman HG. *Coordinated regulation of intracellular pH by two glucose-sensing pathways in yeast.* J Biol Chem. 2018 Feb 16;293(7):2318-2329. doi: [10.1074/jbc.RA117.000422](https://doi.org/10.1074/jbc.RA117.000422). PMID: [29284676](https://pubmed.ncbi.nlm.nih.gov/29284676)
5. Luna LA, Lesecq Z, White KA, Hoang A, Scott DA, Zagnitko O, Bobkov AA, Barber DL, Schiffer JM, Isom DG, Sohl CD. *An acidic residue buried in the dimer interface of isocitrate dehydrogenase 1 (IDH1) helps regulate catalysis and pH sensitivity.* Biochem J. 2020 Aug 28;477(16):2999-3018. doi: [10.1042/BCJ20200311](https://doi.org/10.1042/BCJ20200311). PMID: [32729927](https://pubmed.ncbi.nlm.nih.gov/32729927)
6. Rowe JB, Kapolka NJ, Taghon GJ, Morgan WM, Isom DG. *The evolution and mechanism of GPCR proton sensing.* J Biol Chem. 2021 Jan-Jun;296:100167. doi: [10.1074/jbc.RA120.016352](https://doi.org/10.1074/jbc.RA120.016352). PMID: [33478938](https://pubmed.ncbi.nlm.nih.gov/33478938)
7. Taghon GJ, Rowe JB, Kapolka NJ, Isom DG. *Predictable cholesterol binding sites in GPCRs lack consensus motifs.* Structure. 2021 May 6;29(5):499-506.e3. doi: [10.1016/j.str.2021.01.004](https://doi.org/10.1016/j.str.2021.01.004). PMID: [33508215](https://pubmed.ncbi.nlm.nih.gov/33508215)


