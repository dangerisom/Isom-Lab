---
title: Code
layout: default
permalink: /projects/Superdarks/code/
---

# Project Code Folders

Below are the code sections. Click a folder to view its contents.

{% assign here = page.dir %}

{%- comment -%}
Collect anything (pages or static files) under this directory and then
group by the first path segment after `here` to get *immediate* subfolders.
Works even if a subfolder has no index.md.
{%- endcomment -%}
{% assign under_pages = site.pages | where_exp: "p", "p.dir startswith here and p.dir != here" %}
{% assign under_files = site.static_files | where_exp: "f", "f.path startswith here" %}
{% assign combined = under_pages | concat: under_files %}

{% assign groups = combined
  | group_by_exp: "x", "(x.dir or x.path) | remove_first: here | split: '/' | first" %}

{% assign subfolders = "" | split: "" %}
{% for g in groups %}
  {% assign name = g.name | strip %}
  {% if name != "" %}
    {% assign subfolders = subfolders | push: name %}
  {% endif %}
{% endfor %}
{% assign subfolders = subfolders | uniq | sort %}

{% if subfolders.size == 0 %}
_No subfolders found._
{% else %}
<ul>
  {% for name in subfolders %}
    {% assign subdir = here | append: name | append: "/" %}
    <li>
      📁 <a href="{{ subdir | relative_url }}">{{ name }}/</a>
    </li>
  {% endfor %}
</ul>
{% endif %}
