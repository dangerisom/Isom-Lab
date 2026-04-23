---
layout: default
title: vesicle_triangulator example output — ch1
permalink: /projects/vesicle_triangulator/example_output/ch1/
---
<!-- stager: preserve -->

# vesicle_triangulator — Example output : ch1

Per-channel intermediates for **channel 01** of the example run. These files are the upstream products that feed the joint triangulation in [`../mixed/`](../mixed/); they are kept next to the final outputs so you can inspect how the ch1 vesicle points were detected before being triangulated against ch2.

<ul>
{% assign this_folder = 'github_website/projects/vesicle_triangulator/example_output/ch1/' %}
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

- **`ch1_thresh.png`** — Binary threshold mask for channel 01: the raw pixel-level result of the ch1 threshold setting. This is what you see in the GUI right after clicking "Threshold Ch1".
- **`ch1_components.png`** — Connected-component labeling of `ch1_thresh.png`: each distinct vesicle is colorized with its own component ID so you can spot-check that neighboring vesicles were correctly split (or not merged).
- **`ch1_points.csv`** — Point centroids extracted from the components as `(x, y)` coordinates, one row per vesicle. This is the input to the triangulator — all ch1 vertices in the mesh come from here.
- **`ch1_overlay.png`** — Centroids from `ch1_points.csv` drawn on top of the original microscopy image, so you can verify visually that every detected point lands on a real vesicle.
