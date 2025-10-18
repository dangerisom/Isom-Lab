---
layout: default
title: Superdarks code — folders
---

# 📁 Superdarks: Code folders

<ul>
{%- assign files_root  = "/github_website/projects/Superdarks/code/" -%}
{%- assign pages_root  = "/github_website/projects/Superdarks/code/" -%}
{%- assign subfolder_blob = "" -%}

{%- comment -%} Collect folder names from static files {%- endcomment -%}
{%- for f in site.static_files -%}
  {%- if f.path contains files_root -%}
    {%- assign rel = f.path | remove_first: files_root -%}
    {%- assign parts = rel | split: "/" -%}
    {%- if parts.size > 1 -%}
      {%- assign subfolder_blob = subfolder_blob | append: parts[0] | append: "|" -%}
    {%- endif -%}
  {%- endif -%}
{%- endfor -%}

{%- comment -%} Also collect from pages (to catch folders that have an index.md even without static files) {%- endcomment -%}
{%- for p in site.pages -%}
  {%- if p.path contains pages_root and p.path != page.path -%}
    {%- assign relp = p.path | remove_first: pages_root -%}
    {%- assign partsp = relp | split: "/" -%}
    {%- if partsp.size > 1 -%}
      {%- assign subfolder_blob = subfolder_blob | append: partsp[0] | append: "|" -%}
    {%- endif -%}
  {%- endif -%}
{%- endfor -%}

{%- assign subfolders = subfolder_blob | split: "|" | uniq | sort -%}

{%- if subfolders.size == 0 -%}
  <li><em>No subfolders found.</em></li>
{%- else -%}
  {%- for name in subfolders -%}
    {%- if name != "" -%}
      <li><a href="{{ pages_root | append: name | append: '/' | relative_url }}">{{ name }}</a></li>
    {%- endif -%}
  {%- endfor -%}
{%- endif -%}
</ul>
