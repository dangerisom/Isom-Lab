---
layout: default
title: superdarks code
---

# 🧬 superdarks: Code files : 2-post_query_analysis

<ul>
{% assign this_folder = "/github_website/projects/superdarks/code/2-post_query_analysis/" %}
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

- **`1-upset_analysis_on_filtered_hits.py`** — Builds UpSet intersection diagrams (plus a matching Excel counts sheet) of InterPro-domain membership across the **filtered** hit set, with tqdm progress bars.
- **`2-upset_analysis_on_unfiltered_hits.py`** — Same UpSet analysis run over the **unfiltered** hit set, using `from_indicators` on the boolean membership matrix for robustness at scale.
- **`3-plot_score_contours.py`** — Headless / cluster-friendly plotting pipeline (no CLI, no GUI) that emits 2D and 3D scatter + density-contour plots of coverage vs TM-score vs rank, plus a companion Excel table — optionally per-rank.

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




