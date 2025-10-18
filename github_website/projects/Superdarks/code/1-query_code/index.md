---
layout: default
title: Superdarks code
---

# 🧬 Superdarks: Code files : 1-query_code

<ul>
{% assign this_folder = "/github_website/projects/Superdarks/code/1-query_code/" %}
{% assign github_repo_base = "https://github.com/dangerisom/Isom-Lab/blob/main" %}
{% assign github_tree_base = "https://github.com/dangerisom/Isom-Lab/tree/main" %}

{%- comment -%}
Build two lists from site.static_files:
  • folders: immediate subfolders (e.g., "library")
  • files:   immediate files in this_folder
{%- endcomment -%}
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

<h2>📁 Folders</h2>
<ul>
  {% if folders.size == 0 %}
    <li><em>No subfolders found.</em></li>
  {% else %}
    {% for name in folders %}
      {% if name != "" %}
        <li>
          <!-- Live-site URL (works if the folder has its own index.md) -->
          <a href="{{ this_folder | append: name | append: '/' | relative_url }}">{{ name }}/</a>
          &nbsp;·&nbsp;
          <!-- Always works: GitHub folder view -->
          <a href="{{ github_tree_base }}{{ this_folder }}{{ name }}/" target="_blank" rel="noopener">View on GitHub</a>
        </li>
      {% endif %}
    {% endfor %}
  {% endif %}
</ul>

<h2>📄 Files</h2>
<ul>
  {% if files.size == 0 %}
    <li><em>No files in this folder.</em></li>
  {% else %}
    {% for rel in files %}
      {% if rel != "" %}
        <li>
          📄 <a href="{{ github_repo_base }}{{ this_folder }}{{ rel }}" target="_blank" rel="noopener">{{ rel }}</a>
          – <a href="{{ this_folder | append: rel | relative_url }}" download>Download</a>
        </li>
      {% endif %}
    {% endfor %}
  {% endif %}
</ul>
</ul>

## ⚠️ Disclaimer

This is **research-grade code** provided for academic and exploratory purposes only. It is **not intended for professional or clinical deployment**.



<!--
---
layout: default
title: Superdarks code
---

# 🧬 Superdarks: Code files : 1-query_code

<ul>
{% assign this_folder = "/github_website/projects/Superdarks/code/1-query_code/" %}
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



