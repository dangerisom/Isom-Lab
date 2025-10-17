---
layout: default
title: Superdarks data
---

# 🧬 Superdarks: Data files

<!-- build: {{ site.time }} | rev: {{ site.github.build_revision }} -->

<ul>
{% assign this_folder = "/github_website/projects/Superdarks/data/" %}
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
