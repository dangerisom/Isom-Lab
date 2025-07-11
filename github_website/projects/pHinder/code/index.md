---
layout: default
title: code
---

<ul>
{% assign this_folder = 'github_website/projects/pHinder/code/' %}
{% assign seen_folders = '' %}
{% for file in site.static_files %}
  {% assign path = file.path | remove_first: '/' %}
  {% if path contains this_folder %}
    {% assign rel = path | remove_first: this_folder %}
    {% if rel contains '/' %}
      {% assign subfolder = rel | split: '/' | first %}
      {% unless seen_folders contains subfolder %}
        {% capture seen_folders %}{{ seen_folders }}|{{ subfolder }}{% endcapture %}
        <li><a href="{{ subfolder }}/">{{ subfolder }}/</a></li>
      {% endunless %}
    {% endif %}
  {% endif %}
{% endfor %}
</ul>
