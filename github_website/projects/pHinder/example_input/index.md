---
layout: default
title: pHinder example_input
---
<!-- stager: preserve -->

# 🧬 pHinder — Example input

Example PDB structures used as inputs for the pHinder computational-geometry pipeline. These are the kind of files pHinder consumes directly: standard Protein Data Bank coordinate files that are parsed into `pdbFile` / `pdbFile_cif` objects and then fed to the side-chain topology and convex-hull / Delaunay routines in the `code/` folder.

<ul>
{% assign this_folder = 'github_website/projects/pHinder/example_input/' %}
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

## 📜 Files in this folder

- **`3p0g.pdb`** — Crystal structure of the β2-adrenergic receptor (PDB 3P0G). A 7-transmembrane GPCR used as a canonical example input for pHinder's side-chain topology and buried-ionizable-network analyses.
