---
layout: default
title: superdarks code
---
<!-- stager: preserve -->

# 🧬 superdarks: Code files : 1-query_code/library

<ul>
{% assign this_folder = "/github_website/projects/superdarks/code/1-query_code/library/" %}
{% assign github_repo_base = "https://github.com/dangerisom/Isom-Lab/blob/main" %}

{%- comment -%} Collect immediate subfolders and files {%- endcomment -%}
{% assign folder_blob = "" %}
{% assign file_blob = "" %}

{% for file in site.static_files %}
  {% if file.path contains this_folder %}
    {% assign rel = file.path | remove_first: this_folder %}
    {% assign parts = rel | split: "/" %}
    {% if parts.size > 1 %}
      {% assign folder_blob = folder_blob | append: parts[0] | append: "|" %}
    {% else %}
      {% unless rel == "" %}
        {% assign file_blob = file_blob | append: rel | append: "|" %}
      {% endunless %}
    {% endif %}
  {% endif %}
{% endfor %}

{% assign folders = folder_blob | split: "|" | uniq | sort %}
{% assign files = file_blob | split: "|" | uniq | sort %}

{%- comment -%} List subfolders (folders first) {%- endcomment -%}
{% for name in folders %}
  {% if name != "" %}
    <li>📁 <a href="{{ this_folder | append: name | append: '/' | relative_url }}">{{ name }}/</a></li>
  {% endif %}
{% endfor %}

{%- comment -%} Then list files in this folder {%- endcomment -%}
{% for rel in files %}
  {% if rel != "" %}
    <li>
      📄 <a href="{{ github_repo_base }}{{ this_folder }}{{ rel }}" target="_blank" rel="noopener">{{ rel }}</a>
      – <a href="{{ this_folder | append: rel | relative_url }}" download>Download</a>
    </li>
  {% endif %}
{% endfor %}
</ul>

## 📜 Scripts in this folder

These are the shared computational-geometry / PDB-I/O modules imported by the stage-1 query scripts. They are the Isom-Lab library core and are also used (copied) by other stages and projects.

- **`compGeometry.py`** — Core computational-geometry primitives: tolerance estimation (`geom_tol`), `Vertex` / `Vertex4D` types, `Edge` / `Triangle` simplex classes, plane coefficients (`planeCoefficients3D` / `4D`), general-position tests (`gp2D` / `3D` / `4D`), distance and centroid helpers.
- **`determinants.py`** — Explicit inline 2×2 / 3×3 / 4×4 determinants (no matrix-library dependency), tuned for tight orientation-test loops.
- **`convexHull3D_2_0.py`** — Incremental 3D convex-hull builder plus a `triangulation2D` helper that extracts and labels Delaunay edges for downstream network construction.
- **`convexHull4D_2_22.py`** — Incremental 4D convex-hull builder that uses the horizon-ridge trick, `Simplex1..3` classes, and general-position jostling — the basis for reduced Delaunay-like representations of aligned structures.
- **`pdbFile.py`** — Hand-rolled PDB parser: `PDBfile` (atom / hetatm dicts, chain selection from `-chains.` filenames, gzip-aware I/O), `PseudoAtom`, and the residue-class constants (`NONPOLAR` / `POLAR` / `IONIZABLE` / `ALL_SIDECHAINS` / `ACTIVE`).
- **`pdbFile_cif.py`** — Sibling of `pdbFile.py` for mmCIF input: streams `_atom_site.*` records into the same `PseudoAtom` / residue data model.
- **`pHinderSurface.py`** — pHinder surface reducer: builds ionizable / polar / apolar surface networks from a PDB using `convexHull4D_2_22` + `goFo`, with gzip output and `perf_counter` timers.
- **`goFo.py`** — Recursive "go-forward" walk through a triangulation or network: forbids revisits, carries level / edge bookkeeping, and supports `PSA` filtering and edge-skip lists — the graph-traversal primitive reused across pHinder and superdarks.
- **`ssco.py`** — Secondary-Structure Constrained-Object extractor: `grab_ssco` pulls helical and sheet CA residues from a PDB (via `PseudoAtom` + 4D circumspheres) and optionally writes `_helix`, `_sheet`, or combined PDBs.
- **`3rd_party_requirements.txt`** — External dependencies required by the stage-1 scripts: `bjobs`, `bsub` (LSF), and `numpy`.

## ⚠️ Disclaimer

This is **research-grade code** provided for academic and exploratory purposes only. It is **not intended for professional or clinical deployment**.
