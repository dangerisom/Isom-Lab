---
title: Code
layout: default
permalink: github_website/projects/Superdarks/code/
---

# Project Code Folders

Below are the code sections. Click a folder to view its contents.

{% assign here = page.dir %}
{%- comment -%}
Find all pages underneath this directory, then group by the first path segment
after `code/` to get only immediate subfolders.
{%- endcomment -%}
{% assign children = site.pages
  | where_exp: "p", "p.dir != here and p.dir startswith here"
  | group_by_exp: "p", "p.dir | remove_first: here | split: '/' | first" 
  | sort: "name" %}

{% if children.size == 0 %}
_No subfolders found._
{% else %}
<ul>
{%- for g in children -%}
  {%- assign subdir = here | append: g.name | append: "/" -%}
  <li>
    <a href="{{ subdir | relative_url }}">
      {{ g.name }}
    </a>
  </li>
{%- endfor -%}
</ul>
{% endif %}
