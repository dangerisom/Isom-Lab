---
layout: default
title: superdarks code
---
<!-- stager: preserve -->

# 🧬 superdarks: Code files : 2-post_query_analysis/library

<ul>
{% assign this_folder = "/github_website/projects/superdarks/code/2-post_query_analysis/library/" %}
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

The stage-2 post-query analysis scripts (`1-upset_analysis_on_filtered_hits.py`, `2-upset_analysis_on_unfiltered_hits.py`, `3-plot_score_contours.py`) are self-contained and rely only on mainstream Python packages, so this `library/` folder does not ship bundled modules. For shared Isom-Lab modules (`compGeometry.py`, `pdbFile.py`, `determinants.py`, …) see the sister `library/` folders under `1-query_code/` and `3-post_query_informatics/`.

- **`3rd_party_requirements.txt`** — External dependencies required by the stage-2 scripts: `numpy`, `tqdm`. (UpSet plots additionally require `upsetplot` and `matplotlib`; contour plots require `matplotlib` in a headless backend — pip-install alongside these.)

## ⚠️ Disclaimer

This is **research-grade code** provided for academic and exploratory purposes only. It is **not intended for professional or clinical deployment**.
