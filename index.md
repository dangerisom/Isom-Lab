---
layout: default
title: Isom Lab
---

# Welcome to the Isom Lab on GitHub

We study how cells sense and respond to their environment, with a focus on acidity, dark proteins, and intercellular cooperativity and communication.

## Software

**New to the lab's tools?** Start with **[Setting up your computer](setup/)** — the
one-time, no-typing setup (Mac or Windows) that every tool below assumes.

Our tools are being packaged as installable Python projects. Released tools live in their own
repositories under the [isomlab](https://github.com/isomlab){:target="_blank" rel="noopener"} organization; the remaining
projects are documented here while they are migrated.

### Released

- **[pHinder](https://github.com/isomlab/pHinder){:target="_blank" rel="noopener"}** — ionizable-residue network and surface
  analysis for protein structures · [project page](github_website/projects/pHinder/)
- **[pam_scanning](https://github.com/isomlab/pam_scanning){:target="_blank" rel="noopener"}** — CRISPR/Cas9 guide RNA and
  chimera-insertion primer design across an ORF, with synonymous PAM silencing and BLAST+
  off-target screening · v1.0.1 · Bioconda recipe under review
- **[bioleads](https://github.com/isomlab/bioleads){:target="_blank" rel="noopener"}** — mine biomedical literature for
  enriched terms, co-occurrence networks, and Swanson-style hypothesis leads · v0.1.0 ·
  Bioconda recipe under review
- **[litlog](https://github.com/isomlab/litlog){:target="_blank" rel="noopener"}** — local-first tracker for the papers you
  read and how they interconnect, written as LLM-trainable tagged text · v0.1.1 ·
  Bioconda recipe under review
- **[vesicle_colocalization_quantifier](https://github.com/isomlab/vesicle_colocalization_quantifier){:target="_blank" rel="noopener"}** —
  interactive two-channel vesicle co-localization for microscopy images · v1.0.0 ·
  [project page](github_website/projects/vesicle_colocalization_quantifier/)
- **[bpp_identifier](https://github.com/isomlab/bpp_identifier){:target="_blank" rel="noopener"}** — quantify bridges,
  projections and protrusions at cell–cell boundaries · v1.0.0 ·
  [project page](github_website/projects/bpp_identifier/)
- **[vesicle_triangulator](https://github.com/isomlab/vesicle_triangulator){:target="_blank" rel="noopener"}** — reconstruct the
  topology of vesicle-transfer events by Delaunay triangulation · v1.0.0 ·
  [project page](github_website/projects/vesicle_triangulator/)

### Shared code

- **[isomlab](https://github.com/isomlab/isomlab){:target="_blank" rel="noopener"}** — the three
  modules used by more than one tool: computational geometry and PDB/mmCIF parsing. Each released
  tool ships its own copy, so you do not need this to install or run them; it is the upstream
  source of truth where fixes are made.

### Project pages

- [superdarks](github_website/projects/superdarks/)
