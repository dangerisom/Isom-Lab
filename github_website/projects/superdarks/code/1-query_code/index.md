---
layout: default
title: superdarks code
---
<!-- stager: preserve -->

# 🧬 superdarks: Code files : 1-query_code

<ul>
{% assign this_folder = "/github_website/projects/superdarks/code/1-query_code/" %}
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

- **`1-superDark.py`** — Per-node TM-align runner: aligns the query PDB against a shard of AlphaFold Database predictions on one compute node (driven by `getopt` flags from the job parser).
- **`1-superDark-jobParser.py`** — LSF `bsub` wrapper that distributes `1-superDark.py` across up to 1,000 Triton / Pegasus / BU-SCC nodes, tracks per-node jobs, and supports resubmission of failed shards.
- **`2-collectResults.py`** — Consolidates per-node result objects into a single hit listing, applying configurable `min_tm_score` and `min_residues` thresholds from its CONFIG block.
- **`3-alignResults.py`** — Re-aligns top-ranked hits to the query to produce aligned PDBs suitable for figure generation.
- **`3-alignResults-jobParser.py`** — LSF `bsub` wrapper for `3-alignResults.py` that fans the re-alignment step out across cluster nodes.
- **`4-parseHits.py`** — Parses `rank--coverage-tm--UNIPROT-…-full.pdb` filenames into a structured hit listing (rank, coverage, TM-score, UniProt ID).
- **`4-parseHits-jobParser.py`** — LSF `bsub` wrapper for `4-parseHits.py`.

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




