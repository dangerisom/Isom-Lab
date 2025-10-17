---
title: 1-query_code
layout: default
permalink: /projects/Superdarks/code/1-query_code/
---

# 1-query_code

Browse the Python scripts and subfolders in this section.

{% assign here = page.dir %}

## Subfolders
{%- comment -%}
Find first-level subfolders by looking under this directory
and grouping by the first path segment after `here`.
Works whether the subfolder has its own index page or only static files.
{%- endcomment -%}
{% assign under_pages = site.pages | where_exp: "p", "p.dir startswith here and p.dir != here" %}
{% assign under_files = site.static_files | where_exp: "f", "f.path startswith here" %}
{% assign combined = under_pages | concat: under_files %}
{% assign groups = combined
  | group_by_exp: "x", "(x.dir or x.path) | remove_first: here | split: '/' | first" %}
{% assign subfolders = "" | split: "" %}
{% for g in groups %}
  {% assign name = g.name %}
  {% if name != "" %}
    {% assign subfolders = subfolders | push: name %}
  {% endif %}
{% endfor %}
{% assign subfolders = subfolders | uniq | sort %}

{% if subfolders.size == 0 %}
_No subfolders._
{% else %}
<ul>
  {% for name in subfolders %}
    <li>📁 <a href="{{ here | append: name | append: '/' | relative_url }}">{{ name }}/</a></li>
  {% endfor %}
</ul>
{% endif %}

---

## Python files (immediate)
{%- comment -%}
List only .py files that are directly in this folder (exclude subfolders).
{%- endcomment -%}
{% assign py_all = site.static_files
  | where_exp: "f", "f.path startswith here"
  | where_exp: "f", "f.extname == '.py'" %}

{% assign any_files = false %}
| File | View | Download |
|---|---:|---:|
{% for f in py_all %}
  {% assign rest = f.path | remove_first: here %}
  {% if rest contains '/' %}
    {# skip nested files #}
  {% else %}
| `{{ rest }}` | [View]({{ f.path | relative_url }}) | <a href="{{ f.path | relative_url }}" download>Download</a> |
    {% assign any_files = true %}
  {% endif %}
{% endfor %}

{% unless any_files %}
_No Python files in this folder._
{% endunless %}
