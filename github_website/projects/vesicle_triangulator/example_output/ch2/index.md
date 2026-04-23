---
layout: default
title: vesicle_triangulator example output — ch2
---
<!-- stager: preserve -->

# vesicle_triangulator — Example output : ch2

Per-channel intermediates for **channel 02** of the example run. Same set of artifacts as [`../ch1/`](../ch1/), generated independently by the ch2 threshold / component / point-extraction pass. These ch2 points are combined with the ch1 points in [`../mixed/`](../mixed/) to build the joint triangulation.

<ul>
{% assign this_folder = 'github_website/projects/vesicle_triangulator/example_output/ch2/' %}
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

- **`ch2_thresh.png`** — Binary threshold mask for channel 02: the raw pixel-level result of the ch2 threshold setting.
- **`ch2_components.png`** — Connected-component labeling of `ch2_thresh.png`, with each vesicle colorized by its own component ID.
- **`ch2_points.csv`** — Point centroids extracted from the ch2 components as `(x, y)` coordinates, one row per vesicle. These are the ch2 vertices that enter the triangulation.
- **`ch2_overlay.png`** — Centroids from `ch2_points.csv` drawn on top of the original microscopy image for QA of the ch2 detections.
