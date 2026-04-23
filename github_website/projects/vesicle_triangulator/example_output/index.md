---
layout: default
title: vesicle_triangulator example output
---
<!-- stager: preserve -->

# vesicle_triangulator — Example output

This folder is the full output bundle produced by running **vesicle_triangulator** on the two-channel inputs in [`../input/`](../input/). The tool organizes its outputs into three subfolders plus a top-level parameter log:

1. **`ch1/`** — per-channel-01 artifacts: the thresholded mask, the labeled connected-component image, the point centroids (as CSV), and the centroid overlay.
2. **`ch2/`** — identical set of artifacts for channel-02.
3. **`mixed/`** — the joint outputs that use points from both channels: the full Delaunay triangulation, a trimmed version, edge-colored variants that highlight "mixed" (cross-channel) vs "same" (within-channel) edges, per-vertex mixing histograms, and an Excel summary of the mixing statistics.
4. **`parameter_log.xlsx`** — a single workbook capturing every parameter used for this run (thresholds, colors, trim settings, output paths) so the run is fully reproducible.

<ul>
{% assign this_folder = "/github_website/projects/vesicle_triangulator/example_output/" %}
{% assign github_repo_base = "https://github.com/dangerisom/Isom-Lab/blob/main" %}

{% assign folder_blob = "" %}
{% assign file_blob = "" %}

{% for file in site.static_files %}
  {% if file.path contains this_folder %}
    {% assign rel = file.path | remove_first: this_folder %}
    {% assign parts = rel | split: "/" %}
    {% if parts.size > 1 %}
      {% assign folder_blob = folder_blob | append: parts[0] | append: "|" %}
    {% else %}
      {% unless rel == "" %}
        {% assign file_blob = file_blob | append: rel | append: "|" %}
      {% endunless %}
    {% endif %}
  {% endif %}
{% endfor %}

{% assign folders = folder_blob | split: "|" | uniq | sort %}
{% assign files = file_blob | split: "|" | uniq | sort %}

{% for name in folders %}
  {% if name != "" %}
    <li>📁 <a href="{{ this_folder | append: name | append: '/' | relative_url }}">{{ name }}/</a></li>
  {% endif %}
{% endfor %}

{% for rel in files %}
  {% if rel != "" %}
    <li>
      📄 <a href="{{ github_repo_base }}{{ this_folder }}{{ rel }}" target="_blank" rel="noopener">{{ rel }}</a>
      – <a href="{{ this_folder | append: rel | relative_url }}" download>Download</a>
    </li>
  {% endif %}
{% endfor %}
</ul>

## 📜 Files in this folder

- **`ch1/`** — Channel-01 intermediates: thresholded mask, connected-component labeling, point centroids (CSV), and the centroid overlay. See [`ch1/`](ch1/).
- **`ch2/`** — Channel-02 intermediates: same set of artifacts as `ch1/` but for the second fluorophore. See [`ch2/`](ch2/).
- **`mixed/`** — Joint outputs combining both channels: the full and trimmed Delaunay triangulations, edge-classified overlays ("mixed" cross-channel edges vs "same" within-channel edges), per-vertex mixing histograms, and a mixing-summary workbook. See [`mixed/`](mixed/).
- **`parameter_log.xlsx`** — Per-run parameter log: channel thresholds, colors used to render each channel, triangulation trim settings, and output paths. One workbook per run, used as the audit trail that makes every example reproducible.
