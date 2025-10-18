---
layout: default
title: Superdarks code — 1-query_code
permalink: /github_website/projects/Superdarks/code/1-query_code/
---

# 🧬 Superdarks: Code — 1-query_code

{%- assign ROOT = "/github_website/projects/Superdarks/code/1-query_code/" -%}

<h2>Contents</h2>
<ul>
  {%- comment -%}
  Collect immediate subfolders and immediate files under ROOT.
  A subfolder is detected when a static file path under ROOT has at least one slash after ROOT.
  {%- endcomment -%}

  {%- assign folder_blob = "" -%}
  {%- assign files_blob  = "" -%}

  {%- for f in site.static_files -%}
    {%- if f.path contains ROOT -%}
      {%- assign rel = f.path | remove_first: ROOT -%}
      {%- assign parts = rel | split: "/" -%}
      {%- if parts.size > 1 -%}
        {%- comment -%} It's inside a subfolder; record that subfolder name {%- endcomment -%}
        {%- assign folder_blob = folder_blob | append: parts[0] | append: "|" -%}
      {%- else -%}
        {%- unless rel == "" -%}
          {%- assign files_blob = files_blob | append: rel | append: "|" -%}
        {%- endunless -%}
      {%- endif -%}
    {%- endif -%}
  {%- endfor -%}

  {%- assign folders = folder_blob | split: "|" | uniq | sort -%}
  {%- assign files   = files_blob  | split: "|" | uniq | sort -%}

  {%- comment -%} First: folders (e.g., "library/") {%- endcomment -%}
  {%- for name in folders -%}
    {%- if name != "" -%}
      <li>
        📁 <a href="{{ ROOT | append: name | append: '/' | relative_url }}">{{ name }}/</a>
      </li>
    {%- endif -%}
  {%- endfor -%}

  {%- comment -%} Then: files in this folder {%- endcomment -%}
  {%- for fname in files -%}
    {%- if fname != "" -%}
      <li>
        📄 <a href="{{ (ROOT | append: fname) | relative_url }}" download>{{ fname }}</a>
      </li>
    {%- endif -%}
  {%- endfor -%}

  {%- if folders == empty and files == empty -%}
    <li><em>No items found in this folder.</em></li>
  {%- endif -%}
</ul>

<p style="font-size:0.9em;opacity:0.8">
Tip: To make a folder link (e.g., <code>library/</code>) open a page instead of 404, add an
<code>index.md</code> in that folder with <code>permalink: {{ ROOT }}library/</code>.
</p>



<!--
---
layout: default
title: Superdarks code
---

# 🧬 Superdarks: Code files : 1-query_code

<ul>
{% assign this_folder = "/github_website/projects/Superdarks/code/1-query_code/" %}
{% assign github_repo_base = "https://github.com/dangerisom/Isom-Lab/blob/main" %}

{% for file in site.static_files %}
  {% if file.path contains this_folder %}
    {% assign rel = file.path | remove_first: this_folder %}
    {% unless rel contains "/" %}
      <li>
        📄 <a href="{{ github_repo_base }}{{ file.path }}" target="_blank" rel="noopener">{{ rel }}</a>
        – <a href="{{ file.path | relative_url }}" download>Download</a>
      </li>
    {% endunless %}
  {% endif %}
{% endfor %}
</ul>

## ⚠️ Disclaimer

This is **research-grade code** provided for academic and exploratory purposes only. It is **not intended for professional or clinical deployment**.

-->


