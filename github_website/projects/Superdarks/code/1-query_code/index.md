---
title: Query Code (Python)
layout: default
permalink: /projects/Superdarks/code/1-query_code/
---

# Python files in `1-query_code`

This page lists all `.py` files in this folder for quick viewing or download.

{% assign here = page.dir %}
{% comment %} Collect all static files under this directory that end with .py {% endcomment %}
{% assign py_files = site.static_files
  | where_exp: "f", "f.path startswith here"
  | where_exp: "f", "f.extname == '.py'"
  | sort: "name" %}

{% if py_files.size == 0 %}
_No Python files found._
{% else %}
| File | View | Download |
|---|---:|---:|
{% for f in py_files %}
| `{{ f.name }}` | [View]({{ f.path | relative_url }}) | <a href="{{ f.path | relative_url }}" download>Download</a> |
{% endfor %}
{% endif %}

> Tip: “View” opens the script in your browser. “Download” saves the file directly.
