---
title: Code
layout: default
permalink: /projects/Superdarks/code/
---

# Project Code Folders

Below are the code sections. Click a folder to view its contents.

{% assign here = page.dir %}

{%- comment -%}
Gather pages under this directory (but not the directory itself).
Use 'contains' instead of 'startswith' for GitHub Pages compatibility.
{%- endcomment -%}
{% assign under_pages = site.pages | where_exp: "p", "p.dir != here and p.dir contains here" %}

{%- comment -%}
Gather static files whose *directory* is under this directory.
We derive each file's directory with (f.path | remove: f.name).
{%- endcomment -%}
{% assign under_files = site.static_files 
  | where_exp: "f", "(f.path | remove: f.name) != here and (f.path | remove: f.name) contains here" %}

{%- comment -%}
From pages: take the first segment after 'here'.
{%- endcomment -%}
{% assign page_groups = under_pages 
  | group_by_exp: "p", "p.dir | remove_first: here | split: '/' | first" %}
{% assign names_from_pages = page_groups | map: "name" %}

{%- comment -%}
From files: use each file's directory (path minus filename), then take the first segment after 'here'.
{%- endcomment -%}
{% assign file_groups = under_files 
  | group_by_exp: "f", "(f.path | remove: f.name) | remove_first: here | split: '/' | first" %}
{% assign names_from_files = file_groups | map: "name" %}

{%- comment -%}
Combine, remove blanks, uniq, and sort.
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
_No subfolders found._
{% else %}
<ul>
  {% for name in subfolders %}
    {% assign subdir = here | append: name | append: "/" %}
    <li>📁 <a href="{{ subdir | relative_url }}">{{ name }}/</a></li>
  {% endfor %}
</ul>
{% endif %}
