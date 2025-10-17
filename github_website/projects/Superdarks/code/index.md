---
title: Code
layout: default
permalink: /projects/Superdarks/code/
# IMPORTANT: If your repo keeps these under a folder named "github_website",
# change root below to /github_website/projects/Superdarks/code/
root: /projects/Superdarks/code/
---

# Project Code

Below are the immediate subfolders and files under `{{ page.root }}`.

{% assign root = page.root %}
{% assign files_under = site.static_files | where_exp: "f", "f.path contains root" %}
{% assign pages_under = site.pages | where_exp: "p", "p.path != page.path and p.path contains root" %}

{%- comment -%} Collect immediate subfolder names {%- endcomment -%}
{% assign subfolder_blob = "" %}
{% for f in files_under %}
  {% assign rel = f.path | remove_first: root %}
  {% assign parts = rel | split: "/" %}
  {% if parts.size > 1 %}
    {% assign subfolder_blob = subfolder_blob | append: parts[0] | append: "|" %}
  {% endif %}
{% endfor %}
{% for p in pages_under %}
  {% assign relp = p.path | remove_first: root %}
  {% assign partsp = relp | split: "/" %}
  {% if partsp.size > 1 %}
    {% assign subfolder_blob = subfolder_blob | append: partsp[0] | append: "|" %}
  {% endif %}
{% endfor %}
{% assign subfolders = subfolder_blob | split: "|" | uniq | sort %}

{%- comment -%} Collect immediate top-level files {%- endcomment -%}
{% assign topfiles_blob = "" %}
{% for f in files_under %}
  {% assign rel = f.path | remove_first: root %}
  {% unless rel == "" %}
    {% assign parts = rel | split: "/" %}
    {% if parts.size == 1 %}
      {% assign topfiles_blob = topfiles_blob | append: rel | append: "|" %}
    {% endif %}
  {% endunless %}
{% endfor %}
{% assign topfiles = topfiles_blob | split: "|" | uniq | sort %}

## Subfolders
{% if subfolders.size == 0 %}
_No subfolders found._
{% else %}
<ul>
  {% for name in subfolders %}
    {% if name != "" %}
      <li>
        <a href="{{ root | append: name | append: '/' | relative_url }}">{{ name }}</a>
      </li>
    {% endif %}
  {% endfor %}
</ul>
{% endif %}

## Files
{% if topfiles.size == 0 %}
_No files in this folder._
{% else %}
<ul>
  {% for fname in topfiles %}
    {% assign url = root | append: fname | relative_url %}
    <li>
      <a href="{{ url }}">{{ fname }}</a>
      {% assign ext = fname | split: "." | last | downcase %}
      {% if ext == "py" or ext == "pdb" %}
        &nbsp;—&nbsp;<a href="{{ url }}" download>Download</a>
      {% endif %}
    </li>
  {% endfor %}
</ul>
{% endif %}
