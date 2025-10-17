---
title: 1-query_code
layout: default
permalink: /projects/Superdarks/code/1-query_code/
---

# 1-query_code

Browse the Python scripts and subfolders in this section.

## Subfolders
{% assign here = page.dir %}

{%- comment -%}
Collect first-level subfolder names by looking at both pages and static files
that live under this directory, then grouping by the first path segment after `here`.
{%- endcomment -%}
{% assign under_pages = site.pages | where_exp: "p", "p.dir startswith here and p.dir != here" %}
{% assign under_files = site.static_files | where_exp: "f", "f.path startswith here" %}
{% assign combined = under_pages | concat: under_files %}
{% assign groups = combined
  | group_by_exp: "x", "(x.dir or x.path) | remove_first: here | split: '/' | first" %}
{% assign subfolders = "" | split: "" %}
{% for g in groups %}
  {%- assign name = g.name -%}
  {%- if name != "" -%}
    {%- assign subfolders = subfolders | push: name -%}
  {%- endif -%}
{% endfor %}
{% assign subfolders = subfolders | uniq | sort %}

{% if subfolders.size == 0 %}
_No subfolders._
{% else %}
<ul>
  {% for name in subfolders %}
    <li><a href="{{ here | append: name | append: '/' | relative_url }}">{{ name }}</a></li>
  {% endfor %}
</ul>
{% endif %}

---

## Python files (immediate)
| File | View | Download |
|---|---:|---:|
{% assign py_all = site.static_files
  | where_exp: "f", "f.path startswith here"
  | where_exp: "f", "f.extname == '.py'" %}

{%- comment -%}
Show only files directly in this folder (no slash after removing `here`).
{%- endcomment -%}
{% assign any_files = false %}
{% for f in py_all %}
  {% assign rest = f.path | remove_first: here %}
  {% if rest contains '/' %}
    {# skip nested files #}
  {% else %}
| `{{ f.name }}` | [View]({{ f.path | relative_url }}) | <a href="{{ f.path | relative_url }}" download>Download</a> |
    {% assign any_files = true %}
  {% endif %}
{% endfor %}

{% unless any_files %}
_No Python files in this folder._
{% endunless %}

> “View” opens the script in your browser; “Download” saves it directly.
