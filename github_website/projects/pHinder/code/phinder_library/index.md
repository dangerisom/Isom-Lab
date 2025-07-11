---
layout: default
title: phinder_library_dude
---

<ul>
{% assign this_folder = 'github_website/projects/pHinder/code/phinder_library' %}
{% for file in site.static_files %}
  {% assign path = file.path | remove_first: '/' %}
  {% if path contains this_folder %}
    {% assign rel_path = path | remove_first: this_folder | remove_first: '/' %}
    <li>
      <a href="https://github.com/dangerisom/Isom-Lab/blob/main/{{ path }}">{{ rel_path }}</a>
      — <a href="{{ path | relative_url }}" download>Download</a>
    </li>
  {% endif %}
{% endfor %}
</ul>
