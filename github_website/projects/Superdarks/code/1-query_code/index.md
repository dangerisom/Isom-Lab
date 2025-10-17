---
layout: default
title: 1-query_code
permalink: /projects/Superdarks/code/1-query_code/
---

<ul>
{% assign here = page.dir %} {# e.g., /github_website/projects/Superdarks/code/1-query_code/ #}

{% assign files = site.static_files
  | where_exp: "f", "f.path startswith here"
  | where_exp: "f", "f.extname == '.py'" %}

{% for f in files %}
  {% assign rest = f.path | remove_first: here %}
  {% unless rest contains '/' %}
    <li>
      📄 <a href="{{ f.path | relative_url }}" target="_blank">{{ rest }}</a>
      – <a href="{{ f.path | relative_url }}" download>Download</a>
    </li>
  {% endunless %}
{% endfor %}
</ul>
