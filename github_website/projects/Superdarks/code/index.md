---
title: Code
layout: default
permalink: /projects/Superdarks/code/
---

# Project Code Folders

Below are the code sections. Click a folder to view its contents.

{%- comment -%}
Use page.url as the canonical base. Emulate "startsWith" via slice compare.
This lists *immediate* subfolders under /projects/Superdarks/code/.
Works even if a subfolder has no index.md.
{%- endcomment -%}
{% assign base = page.url %}

{%- comment -%}
Collect pages under base (exclude the base itself)
{%- endcomment -%}
{% assign under_pages = site.pages
  | where_exp: "p", "(p.url | slice: 0, base.size) == base and p.url != base" %}

{%- comment -%}
Collect static files under base (exclude the base itself)
For each file, work with its directory path = f.path without filename.
{%- endcomment -%}
{% assign under_files = site.static_files
  | where_exp: "f", "((f.path | slice: 0, base.size) == base) and f.path != base" %}

{%- comment -%}
Get first segment after base for pages (using url)
{%- endcomment -%}
{% assign page_groups = under_pages
  | group_by_exp: "p", "p.url | remove_first: base | split: '/' | first" %}
{% assign names_from_pages = page_groups | map: "name" %}

{%- comment -%}
Get first segment after base for files (using directory path)
{%- endcomment -%}
{% assign file_groups = under_files
  | group_by_exp: "f", "(f.path | remove: f.name) | remove_first: base | split: '/' | first" %}
{% assign names_from_files = file_groups | map: "name" %}

{%- comment -%}
Combine, clean, uniq, sort
{%- endcomment -%}
{% assign all_names = names_from_pages | concat: names_from_files %}
{% assign cleaned = "" | split: "" %}
{% for n in all_names %}
  {% assign name = n | to_s | strip %}
  {% if name != "" %}
    {% assign cleaned = cleaned | push: name %}
  {% endif %}
{% endfor %}
{% assign subfolders = cleaned | uniq | sort %}

{% if subfolders.size == 0 %}
_No subfolders found under {{ base }}._
{% else %}
<ul>
  {% for name in subfolders %}
    {% assign subdir = base | append: name | append: "/" %}
    <li>📁 <a href="{{ subdir | relative_url }}">{{ name }}/</a></li>
  {% endfor %}
</ul>
{% endif %}
