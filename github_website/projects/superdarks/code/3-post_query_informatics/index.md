---
layout: default
title: superdarks code
---

# 🧬 superdarks: Code files : 3-post_query_informatics

<ul>
{% assign this_folder = "/github_website/projects/superdarks/code/3-post_query_informatics/" %}
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

- **`1-listing_parse.py`** — Parses `rank / coverage / TM-score / UniProt-ID` tuples out of the `rank--coverage-tm--UNIPROT-…-full.pdb` filenames produced upstream.
- **`2-uniprot_parse.py`** — Streaming parser for Swiss-Prot / TrEMBL `.dat` flat-file entries, one protein record at a time (memory-safe).
- **`3-listing_uniprot_merge.py`** — Streams the full Swiss-Prot and TrEMBL `.dat.gz` archives via `pigz` / `gzip` and merges them with the hit listing into a unified TSV (gene, organism, superkingdom, and protein-function columns).
- **`4-uniprot_interpro_merge.py`** — Joins the merged UniProt TSV with InterPro's `protein2ipr.dat.gz`, emitting one combined TSV per hit (UPPERCASE-normalized for downstream matching).
- **`5-gpcr_filtering.py`** — Pandas-based 7TM + GPCR keyword filter over the InterPro-enriched TSV (with a special `GCR-\d+` family rule), writing per-match Excel workbooks up to the 1,048,576-row sheet limit.
- **`6-sequence_chunkify_and_collect.py`** — Rewrites FASTA headers into the canonical `superkingdom|uniprot_id|gene|organism|rank` form and chunks sequences into downstream-ready files.
- **`7-clusterupdate_branch/`** — Optional pre-BLASTP reduction: splits the FASTA by superkingdom and runs MMseqs2 `linclust` + `clusterupdate` to produce priority representatives (Archaea → Bacteria → Eukaryota).
- **`8-chunkify_and_blastp.py`** — Splits one FASTA into N whole-record chunks and submits one all-vs-all BLASTP LSF job per chunk (no CLI — edit the CONFIG block).
- **`9-blastp_chunks_to_network_fast.py`** — Unified network builder (FAST in-process, or TURBO / LSF cluster) that turns per-chunk BLASTP TSVs into gzipped `nodes` + `edges` tables ready for Cytoscape.
- **`10-goFo_one_node.py`** — Isom-Lab "go-forward" BFS walk over the thresholded BLAST similarity graph (by % identity plus e-value or bitscore), producing a no-rings, strongest-edges-first subgraph.
- **`11-reduce_for_cytoscape_strict.py`** — Reduces a very large network to a Cytoscape-friendly subset by keeping the largest components, extracting a Maximum Spanning Forest backbone, and pruning by node-importance score.
- **`12-clustalo_pairwise_traceback_map_con_to_structure.py`** — Clustal-Omega pairwise traceback that maps sequence conservation onto a reference PDB, enabling 3D visualization of evolutionary signal.

## ⚠️ Disclaimer

This is **research-grade code** provided for academic and exploratory purposes only. It is **not intended for professional or clinical deployment**.


<!--
---
layout: default
title: superdarks code
---

# 🧬 superdarks: Code files : 1-query_code

<ul>
{% assign this_folder = "/github_website/projects/superdarks/code/1-query_code/" %}
{% assign github_repo_base = "https://github.com/dangerisom/Isom-Lab/blob/main" %}

{% for file in site.static_files %}
  {% if file.path contains this_folder %}
    {% assign rel = file.path | remove_first: this_folder %}
    {% unless rel contains "/" %}
      <li>
        📄 <a href="{{ github_repo_base }}{{ file.path }}" target="_blank" rel="noopener">{{ rel }}</a>
        – <a href="{{ file.path | relative_url }}" download>Download</a>
      </li>
    {% endunless %}
  {% endif %}
{% endfor %}
</ul>

## ⚠️ Disclaimer

This is **research-grade code** provided for academic and exploratory purposes only. It is **not intended for professional or clinical deployment**.
-->




