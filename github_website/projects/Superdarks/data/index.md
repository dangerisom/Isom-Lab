---
title: Data
layout: default
permalink: /projects/Superdarks/data/
# Pretty page URL:
root_pages: /projects/Superdarks/data/
# Real static file location (note github_website):
root_files: /github_website/projects/Superdarks/data/
---

{% assign P = page.root_pages %}
{% assign F = page.root_files %}

{%- comment -%} List files from F, link pages/subfolders from P {%- endcomment -%}
{% assign files_under = site.static_files | where_exp: "f", "f.path contains F" %}
{% assign pages_under = site.pages | where_exp: "p", "p.path != page.path and p.path contains P" %}

<!-- build subfolders same as before, but link with P -->
<!-- build top-level files same as before, but link with F -->

<!-- Example of a file link -->
{% for f in files_under %}
  {% assign rel = f.path | remove_first: F %}
  {% assign parts = rel | split: "/" %}
  {% if parts.size == 1 %}
  {% assign file_url = F | append: rel | relative_url %}
  <li>
    <a href="{{ file_url }}">{{ rel }}</a>
    {% assign ext = rel | split: "." | last | downcase %}
    {% if ext == "py" or ext == "pdb" %}
      — <a href="{{ file_url }}" download>Download</a>
    {% endif %}
  </li>
  {% endif %}
{% endfor %}
