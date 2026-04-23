# Isom Lab — github_website

Source tree for the public Isom Lab website at **https://dangerisom.github.io/Isom-Lab/**.
Jekyll renders the top-level `Isom-Lab/` repo with the `minima` theme; this
`github_website/` folder is where every published *project* lives. A project is a
self-contained scientific tool (one entry script or a staged pipeline) plus the data
and documentation that ships with it.

The layout is fixed by `stage_project_for_website.py` in the private tooling repo
(`IsomLabPrivate/pythonLibrary/`), and every page here is either hand-written Markdown
or auto-generated from a Python script's docstring, argparse definition, and
dependency list.

All code in this tree is licensed under **MPL 2.0**; the header is stamped into every
`.py` file automatically by the stager.

---

## Site structure

```
github_website/
├── README.md                              (this file)
├── initial_prompt.txt                     (2025 design notes for the Jekyll site)
└── projects/
    ├── Superdarks/
    ├── bpp_identifier/
    ├── pHinder/
    ├── vesicle_colocalization_quantifier/
    └── vesicle_triangulator/
```

Each `projects/<name>/` folder follows the same skeleton:

```
<name>/
├── index.md        ← Jekyll landing page (front-matter + overview)
├── README.md       ← (optional) plain-Markdown twin of index.md
├── code/
│   ├── index.md    ← one row per file under code/; hand-written captions OK
│   ├── <file>.py   ← entry script(s) and library modules
│   └── <sub>/      ← optional subfolders with their own index.md
└── data/
    └── index.md    ← one row per file under data/; hand-written captions OK
```

The stager regenerates `index.md` and `README.md` at the project root from each
entry script's docstring on every run — **unless** the existing file contains the
sentinel comment `<!-- stager: preserve -->`, in which case the stager leaves it
alone. Drop that HTML comment just below the YAML front matter on any hand-edited
landing page you want to protect. The `index.md` files under `code/` and `data/`
are **never overwritten** regardless — the stager only creates them when missing,
so hand-written file captions survive every rerun.

---

## Project overviews

### `pHinder/` — pH-dependent interaction networks and side-chain topology

The flagship Isom-Lab tool. **pHinder** classifies every side chain in a protein as
**core**, **margin**, or **exposed** relative to a Cα-based Delaunay surface, detects
contiguous networks of side chains (especially ionizable ones), and writes PDB + CONECT
output files that visualize each network in ChimeraX / PyMOL.

Two entry points ship side by side:

- `phinder_main_gui.py` — full Tkinter interface (custom `DynamicOptionWidget`,
  `AminoAcidSelectionWidget`, `FilePathWidget`, collapsible panes, redirected stdout
  piped into a live `TerminalOutputWidget`).
- `phinder_command_line.py` — argparse CLI exposing the same parameters (residue set,
  topology vs surface vs side-chain classification, interface detection, core / margin
  cutoffs, max edge length, virtual-screening options, reduced network representation).

The algorithm itself lives in `pHinder_7_0.py`, which composes the geometry primitives
in `compGeometry.py` and `determinants.py`, the 3D and 4D hull engines in
`convexHull3D_1_0.py` and `convexHull4D_2_22.py`, the molecular-surface code in
`pHinderSurface.py`, the goFo recursive traversal in `goFo.py`, the network minimizer
in `minimizeNetworks.py`, the PDB + mmCIF parser in `pdbFile.py` / `pdbFile_cif.py`,
sphere sampling in `sphere.py`, and PDB / CONECT writers in `writeFunctions.py`. The
rest of the folder (`collapsible_toggle.py`, `dynamic_option_widget*.py`,
`file_open.py`, `input_widgets.py`, `terminal_output_widget.py`) is Tkinter widget
scaffolding used only by the GUI.

Beyond the core topology calculation, pHinder also supports protein–protein interface
detection and ligand-binding-cavity prediction. The method has been published in *Mol
Cell* (2013) and *PNAS* (2015); the project's `index.md` lists the full citation set.

### `Superdarks/` — large-scale structural alignment against the AlphaFold Database

A four-stage pipeline for finding structural homologs ("darks", hence *Superdarks*) in
the **AlphaFold Database** (≈214 M predictions) starting from one or more query PDBs.
It distributes TM-align jobs across 1000 HPC nodes on University-of-Miami Triton (or
Pegasus / BU SCC), then post-processes the hits through progressively richer
annotation and network-building stages. Each subfolder is an ordered pipeline step:

**`1-query_code/`** drives the structural-alignment sweep. `1-superDark.py` is the
per-node TM-align runner; `1-superDark-jobParser.py` is its job-submission wrapper.
`2-collectResults.py` consolidates per-node output into a single result object.
`3-alignResults.py` (+ its `-jobParser` companion) re-aligns top hits for figure
generation. `4-parseHits.py` (+ `-jobParser`) walks the aligned-hits folder and
extracts score-tagged PDB files using coverage × TM-score thresholds.

**`2-post_query_analysis/`** summarizes hits across domains/queries.
`1-upset_analysis_on_filtered_hits.py` and `2-upset_analysis_on_unfiltered_hits.py`
build UpSet-diagram intersection counts (via `upsetplot`) over domain membership.
`3-plot_score_contours.py` is a headless plotting pipeline (no CLI, no GUI) that emits
2D/3D scatter and density-contour plots of coverage vs TM vs rank plus an Excel table.

**`3-post_query_informatics/`** enriches hits with annotations and builds an
interaction network.
- `1-listing_parse.py` / `2-uniprot_parse.py` — parse `ls -l` listings and Swiss-Prot
  / TrEMBL `.dat` entries into (rank, coverage, TM) and (AC, gene, organism, …) dicts.
- `3-listing_uniprot_merge.py` — streams both Swiss-Prot and TrEMBL through
  `pigz` / `gzip`, merges against a hit listing, emits a unified TSV.
- `4-uniprot_interpro_merge.py` — joins the UniProt TSV with InterPro annotations in
  UPPERCASE.
- `5-gpcr_filtering.py` — pandas / openpyxl TSV filter with 7TM+GPCR keyword matching
  and a GCR-family rule.
- `6-sequence_chunkify_and_collect.py` — filters a FASTA by a UniProt-ID list and
  writes overlapping pair-FASTAs, with a UniProt-accession regex fallback.
- `8-chunkify_and_blastp.py` / `9-blastp_chunks_to_network_fast.py` — chunk one FASTA,
  submit BLASTP jobs via LSF `bsub`, then build a nodes/edges network (FAST in-process
  or TURBO/LSF) gzipped for Cytoscape.
- `10-goFo_one_node.py` — thresholded go-forward walk on the BLAST similarity graph
  (filter by `%identity`, `evalue`, or bitscore).
- `11-reduce_for_cytoscape_strict.py` — reduces a very large network to a
  Cytoscape-friendly subset using stock Python 3 (no extra deps).
- `12-clustalo_pairwise_traceback_map_con_to_structure.py` — Clustal-Omega pairwise
  traceback mapping of conservation onto a reference PDB structure.

**`4-python_wrappers/`** holds standalone wrappers used at any stage:
`foldseek_individual_query_gui.py` is a Tkinter front-end for running Foldseek
`createdb` + search with stale-index detection; `fpocket_python_wrapper.py` is a thin
subprocess wrapper around `fpocket` that gathers its `*_out/` directories.

### `bpp_identifier/` — Bridge / Projection / Protrusion image quantifier

Single-file Tkinter app (`bpp_identifier_v22_1.py`) that loads a microscopy overlay
image, applies a user-adjustable threshold / contour pipeline (OpenCV), and quantifies
bridges, projections, and protrusions at a cell boundary. The window bundles debugging
toggles (`enable_interface_islands`, `enable_green_fragment_to_interface`), a radial
search radius for skeleton-midline connection detection, and a per-run save directory
that is time-stamped so repeat runs never collide. Dependencies are `cv2`, `numpy`,
and `Pillow`.

### `vesicle_colocalization_quantifier/` — two-channel vesicle co-localization

Single-file Tkinter app (`vesicle_coloc_overlap_v2_1.py`) that takes two microscopy
channels and measures pixel-level co-localization. The GUI lets the user set
independent low/mid/high thresholds per channel, configure contour area bounds and
line width, and exports counts (`total1`, `total2`, `overlap_count`), fractional
overlaps (`frac1`, `frac2`), and yellow-outlined overlap contours. Dependencies are
`cv2`, `numpy`, `csv`.

### `vesicle_triangulator/` — 3D Delaunay triangulation of vesicle transfer events

Single-file Tkinter app (`vesicle_transfer_triangulation_v10_1.py`) that reads CSV
coordinates plus component images, constructs a 3D convex-hull / Delaunay
triangulation via the Isom-Lab `compGeometry` / `convexHull3D_2_1` modules, and writes
the result as both an overlay image and an auto-fit Excel workbook (column widths
adjusted via `openpyxl.utils.get_column_letter`). Dependencies are `cv2`, `numpy`,
`compGeometry`, `convexHull3D_2_1`, `pandas`, `matplotlib`, `Pillow`, `openpyxl`, and
the pHinder-pipeline geometry library.

---

## Conventions

**Licensing.** Every `.py` file opens with the MPL 2.0 header and an Isom-Lab / Daniel
G. Isom copyright line. The stager in `IsomLabPrivate/pythonLibrary/_collect_common.py`
stamps this automatically; do not remove it.

**Entry-point styles.** The projects mix three paradigms:

- *Tkinter GUIs* (`phinder_main_gui.py`, `bpp_identifier_v22_1.py`,
  `vesicle_coloc_overlap_v2_1.py`, `vesicle_transfer_triangulation_v10_1.py`,
  `foldseek_individual_query_gui.py`) — run with `python <file>.py`.
- *argparse CLIs* (`phinder_command_line.py`) — parse flags, expose `--help`.
- *Configure-and-run scripts* — most of `Superdarks/code/` edits a CONFIG block near
  the top (or takes `getopt` flags for HPC submission) rather than a full CLI. This is
  deliberate: cluster jobs are driven by `bsub` / Slurm wrappers, not interactive
  input.

**Jekyll front matter.** Every `index.md` starts with a YAML block (`layout: default`,
`title: …`, and usually `permalink: …`). Internal links use the
`{{ site.baseurl }}/github_website/projects/<name>/...` pattern so `url` / `baseurl`
changes in `_config.yml` propagate automatically. The landing `index.md` at
`Isom-Lab/index.md` links projects individually — add a new line there to feature a
project.

**Dependencies.** Each project's `code/` has a `3rd_party_requirements.txt`. Core
dependency set across the tree: `numpy`, `cv2` (OpenCV), `pandas`, `openpyxl`,
`matplotlib`, `tqdm`, `Pillow`. Larger pipelines additionally require `scipy`,
`primer3-py`, `rdkit`, `requests` (UniProt), `upsetplot`, and cluster-side binaries
(`TMalign`, `fpocket`, `foldseek`, `blastp`, `pigz`, Clustal-Omega, LSF `bsub`).

**File paths.** Many `Superdarks/` scripts hard-code cluster paths under
`/projectnb/isomlab/dan/…` (BU SCC) or `/nethome/dgi5/…` (Triton). Before running on a
new machine, search the header of the script for the CONFIG block and adjust paths
and `cluster` / `queue` variables.

**Temporary files and `__pycache__`.** The site-root `.gitignore` excludes `.DS_Store`,
`_site/`, `.jekyll-cache/`, `.jekyll-metadata`, `.sass-cache/`, `__pycache__/`, and
`*.pyc`. Stage new projects from a clean tree (the stager skips `__pycache__/` and
`library/` during recursive folder input).

---

## Adding a new project

1. In `IsomLabPrivate/pythonLibrary/`, run either the CLI or GUI stager:
   `python stage_project_for_website.py <entry.py> <library_folder> <project_name>`
   (or the `_gui.py` twin).
2. That writes a fresh `projects/<project_name>/` folder here with `index.md`,
   `README.md`, `code/`, and `data/`, and stamps every `.py` with the MPL header.
3. Hand-edit the per-file captions under `code/index.md` and `data/index.md` — those
   are preserved on subsequent reruns.
4. Add a bullet to `Isom-Lab/index.md` so the new project shows up on the landing
   page.

---

*Maintained by Daniel G. Isom. The per-project `index.md` files at project roots are
auto-regenerated from each entry script's docstring, argparse definition, and
dependency list by `_readme_gen.py`; do not hand-edit them. Long-form narrative
belongs in `code/index.md` or `data/index.md`.*
