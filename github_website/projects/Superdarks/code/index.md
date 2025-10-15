---
layout: default
title: code
---

<ul>
{% assign this_folder = 'github_website/projects/Superdarks/code/' %}
{% assign github_repo_base = 'https://github.com/dangerisom/Isom-Lab/blob/main/' %}

{% assign subdirs = "" | split: "" %}

{% for file in site.static_files %}
  {% assign path = file.path | remove_first: '/' %}
  {% if path contains this_folder %}
    {% assign rel_path = path | remove_first: this_folder %}
    {% assign parts = rel_path | split: '/' %}
    {% if parts.size > 1 %}
      {% assign subdir = parts[0] %}
      {% unless subdirs contains subdir %}
        {% assign subdirs = subdirs | push: subdir %}
      {% endunless %}
    {% endif %}
  {% endif %}
{% endfor %}

{% assign sorted_subdirs = subdirs | sort %}

{% for folder in sorted_subdirs %}
  <li>
    📁 <a href="{{ github_repo_base }}{{ this_folder }}{{ folder }}/" target="_blank">{{ folder }}</a>
  </li>
{% endfor %}
</ul>
