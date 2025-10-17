---
layout: default
title: code
---

<ul>
{% assign this_folder = 'github_website/projects/Superdarks/code/1-query_code/' %}
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


---

# 🧬 Superdarks: Structural Analysis Tools

Welcome to the **Superdarks** codebase. This folder contains Python code used in our lab's research for protein structural analysis. More to come here when the manuscript publishes. 

## ⚠️ Disclaimer

This is **research-grade code** provided for academic and exploratory purposes only. It is **not intended for professional or clinical deployment**.

- The code was developed and extended over many years to support diverse use cases in computational structural biology.
- You may encounter **commented sections**, **experimental blocks**, or **legacy fragments** — these reflect the evolving nature of scientific coding.
- **Use at your own risk.** While we strive for accuracy, this code has not undergone rigorous software engineering review.

