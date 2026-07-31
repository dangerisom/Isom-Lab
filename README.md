# Isom Lab on GitHub

Source for the Isom Lab code website, plus the project code staged alongside it.

**Rendered site: https://dangerisom.github.io/Isom-Lab/**

We study how cells sense and respond to their environment, with a focus on acidity, dark
proteins, and intercellular cooperativity and communication.

## What is in this repository

```
index.md, about.md, publications.md   the website pages (Jekyll / GitHub Pages)
github_website/projects/              per-project pages, code, and example data
```

| Project | Contents here |
|---|---|
| [superdarks](github_website/projects/superdarks/) | code, example input — structural-homology search across the AlphaFold DB, plus the post-query informatics pipeline |
| [bpp_identifier](github_website/projects/bpp_identifier/) | now at [isomlab/bpp_identifier](https://github.com/isomlab/bpp_identifier); staged copy kept here for provenance |
| [vesicle_colocalization_quantifier](github_website/projects/vesicle_colocalization_quantifier/) | now at [isomlab/vesicle_colocalization_quantifier](https://github.com/isomlab/vesicle_colocalization_quantifier); staged copy kept here for provenance |
| [vesicle_triangulator](github_website/projects/vesicle_triangulator/) | now at [isomlab/vesicle_triangulator](https://github.com/isomlab/vesicle_triangulator); staged copy kept here for provenance |
| [pHinder](github_website/projects/pHinder/) | code, example input — the manuscript version, kept for archival completeness; development continues at [isomlab/pHinder](https://github.com/isomlab/pHinder) |

## Released tools

Packaged, installable tools live in their own repositories under the
[**isomlab**](https://github.com/isomlab) organization:

- [**pHinder**](https://github.com/isomlab/pHinder) — ionizable-residue network and surface
  analysis for protein structures
- [**pam_scanning**](https://github.com/isomlab/pam_scanning) — CRISPR/Cas9 guide RNA and
  chimera-insertion primer design across an ORF
- [**bioleads**](https://github.com/isomlab/bioleads) — mine biomedical literature for enriched
  terms, co-occurrence networks, and Swanson-style hypothesis leads
- [**litlog**](https://github.com/isomlab/litlog) — local-first tracker for the papers you read
  and how they interconnect
- [**vesicle_colocalization_quantifier**](https://github.com/isomlab/vesicle_colocalization_quantifier) — interactive two-channel vesicle co-localization
- [**bpp_identifier**](https://github.com/isomlab/bpp_identifier) — bridges, projections and protrusions at cell–cell boundaries
- [**vesicle_triangulator**](https://github.com/isomlab/vesicle_triangulator) — topology of vesicle-transfer events
- [**isomlab**](https://github.com/isomlab/isomlab) — the geometry and PDB/mmCIF modules shared
  by more than one tool. The released tools vendor their own copies, so this is the upstream
  source of truth rather than an install-time dependency.

## Archived release

The state of this repository as cited in the associated manuscript is preserved at the
[**`nature-code-2026`**](https://github.com/dangerisom/Isom-Lab/tree/nature-code-2026) tag,
which contains all five project folders including the self-contained pHinder code.

Archived at Zenodo:

- **Code** — [10.5281/zenodo.21175708](https://doi.org/10.5281/zenodo.21175708) (MIT)
- **Datasets S1–S4** — [10.5281/zenodo.21178101](https://doi.org/10.5281/zenodo.21178101) (CC BY 4.0)

Both are concept DOIs and always resolve to the latest version. The full citation will be
added on publication.

## Lab website

https://www.isomlab.com

## License

MIT © Daniel G. Isom (Isom Lab)
