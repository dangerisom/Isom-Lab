---
layout: default
title: Superdarks data
permalink: /projects/Superdarks/data/
---

<ul>
{% assign this_folder = 'github_website/projects/Superdarks/data/' %}
{% assign github_repo_base = 'https://github.com/dangerisom/Isom-Lab/blob/main/' %}
{% for file in site.static_files %}
  {% assign path = file.path | remove_first: '/' %}
  {% if path contains this_folder %}
    {% assign rel = path | remove_first: this_folder %}
    {% unless rel contains '/' %}
      <li>
        📄 <a href="{{ github_repo_base }}{{ this_folder }}{{ rel }}" target="_blank">{{ rel }}</a>  
        – <a href="{{ site.baseurl }}/{{ path }}" download>Download</a>
      </li>
    {% endunless %}
  {% endif %}
{% endfor %}
</ul>