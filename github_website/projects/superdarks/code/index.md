---
layout: default
title: superdarks code
---

# 🧬 superdarks: Code files

The superdarks codebase is organized as a **four-stage pipeline**. Each subfolder below is one stage; click into a stage to browse its scripts and to read the stage-specific notes.

1. **`1-query_code/`** — per-node TM-align runner and its job-submission wrapper, result consolidation, top-hit re-alignment, and coverage × TM-score hit parsing.
2. **`2-post_query_analysis/`** — UpSet diagram intersection analysis over filtered and unfiltered hits, plus headless 2D/3D score-contour plots.
3. **`3-post_query_informatics/`** — UniProt / TrEMBL / InterPro enrichment, 7TM+GPCR keyword filtering, BLASTP network construction, go-forward graph traversal, and Clustal-Omega conservation-to-structure mapping.
4. **`4-python_wrappers/`** — standalone wrappers used at any stage: a Foldseek GUI front-end and an `fpocket` subprocess wrapper.

<ul>
{% assign this_folder = "/github_website/projects/superdarks/code/" %}
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

## ⚠️ Disclaimer

This is **research-grade code** provided for academic and exploratory purposes only. It is **not intended for professional or clinical deployment**.

- The code was developed and extended over many years to support diverse use cases in computational structural biology.
- You may encounter **commented sections**, **experimental blocks**, or **legacy fragments** — these reflect the evolving nature of scientific coding.
- **Use at your own risk.** While we strive for accuracy, this code has not undergone rigorous software engineering review.

## 📥 Getting Started

Most scripts in `superdarks/code/` are driven either by a **CONFIG block** at the top of the file (edit in place) or by **`getopt` flags** for HPC submission via LSF `bsub`. Before running on a new machine:

1. Open the script's header and set `cluster`, `queue`, and any absolute paths under `/projectnb/isomlab/…` or `/nethome/…` to match your environment.
2. Confirm required cluster-side binaries are on `PATH`: `TMalign`, `blastp`, `clustalo`, `foldseek`, `fpocket`, `pigz`, and `bsub`.
3. Install Python dependencies listed in each stage's `3rd_party_requirements.txt`.

See the project landing page for the full method overview and citations.
