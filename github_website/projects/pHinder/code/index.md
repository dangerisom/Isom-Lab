---
layout: default
title: code
---

<ul>
{% assign this_folder = 'github_website/projects/pHinder/code/' %}
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

# 🧬 pHinder: Structural Analysis Tools

Welcome to the **pHinder** codebase. This folder contains Python code used in our lab's research for protein structural analysis, sidechain topology and spatial classification, surface calculation, and void detection in protein structures for ligand/drug binding.

## ⚠️ Disclaimer

This is **research-grade code** provided for academic and exploratory purposes only. It is **not intended for professional or clinical deployment**.

- The code was developed and extended over many years to support diverse use cases in computational structural biology.
- You may encounter **commented sections**, **experimental blocks**, or **legacy fragments** — these reflect the evolving nature of scientific coding.
- **Use at your own risk.** While we strive for accuracy, this code has not undergone rigorous software engineering review.

## 📥 Getting Started

Once you have downloaded the files in the `code/` folder to a folder on your computer:

1. Ensure you have Python 3.8+ installed.
2. Install any necessary dependencies (see `requirements.txt` if provided).
3. Prepare your input files in `.pdb` or `.cif` format and place them in an accessible directory.

### Command Line Usage

To run pHinder from the command line, use the `phinder_command_line.py` script:

```bash
python phinder_command_line.py path/to/your_structure.pdb --chains A B --topology-calculation --sidechain-classification --interface-classification
```

You can customize your run by passing additional flags. Use `--help` to see all available options:

```bash
python phinder_command_line.py --help
```

### GUI Mode Usage

To run pHinder with a graphical interface, use the `phinder_main_gui.py` script:

```bash
python phinder_main_gui.py
```

This will launch the GUI, allowing you to select input files and parameters through the interface.

## 📄 Citation

If you use **pHinder** in your research or publications, please **cite the relevant work** from our lab. This supports continued development and scientific recognition.

Thank you for your interest in our research tools!