---
layout: default
title: Project bpp_identifier
---
<!-- stager: preserve -->

# Project bpp_identifier

This project quantifies **bridges, projections, and protrusions** at cell–cell and cell–substrate boundaries in microscopy overlay images using interactive contour analysis.

- [code](code/)
- [data](data/)

## 🔬 What bpp_identifier Does

**bpp_identifier** is a **Tkinter-based image-analysis application** that loads a microscopy overlay, applies a user-adjustable threshold / contour pipeline (OpenCV), and automatically detects and quantifies three classes of boundary morphology: **bridges** (narrow links spanning gaps), **projections** (outward extensions from a cell boundary), and **protrusions** (localized bulges or blebs). It is designed for biologists who want reproducible, per-image quantification without writing code — every run is time-stamped into its own save directory so repeat analyses never collide.

## 🔍 Key Capabilities

- **Interactive Image Loading and Display**  
  Loads microscopy overlays (PNG, TIFF, JPEG) via a file-picker and renders them in-app with Pillow and OpenCV; the display scales to fit while preserving aspect ratio.

- **Adjustable Threshold and Contour Pipeline**  
  Exposes **threshold_min**, **threshold_max**, and hole-area parameters as GUI fields, so the user can dial in segmentation quality without re-running the script or editing source.

- **Skeleton-Midline Connection Detection**  
  Uses a configurable radial search radius (`max_r`) to walk the skeleton midline and detect bridging connections between otherwise-separate contour regions.

- **Separate Hole and Non-Hole Paths**  
  Distinguishes between interior holes (e.g., closed bleb interiors) and exterior contours using independent area thresholds (`min_hole_area_holes`, `hard_min_large_green_area`), avoiding the common failure mode where hole artifacts dominate the count.

- **Debugging Toggles**  
  `enable_interface_islands` and `enable_green_fragment_to_interface` are exposed as runtime flags so the user can step through intermediate stages of the pipeline when tuning parameters for a new image set.

- **Edge-Margin Filtering**  
  An `edge_margin` parameter in pixels draws an inner bounding box that ignores contours touching the image border, preventing imaging artifacts from inflating counts.

- **Per-Run Timestamped Save Directories**  
  Every analysis writes into a dated subfolder so the output directory accumulates a clean audit trail of every parameter sweep.

## 🛠️ Additional Functionalities

- **Original and Processed Overlay Display**  
  Side-by-side display of the input image and the annotated output, making it easy to spot-check segmentation quality before exporting counts.

- **Clean Shutdown**  
  Tracks all `after` callbacks and gracefully cancels them on close, so the GUI never leaves Tkinter threads hanging when the user closes the window mid-analysis.

- **Portable, Minimal Dependencies**  
  Requires only `cv2`, `numpy`, and `Pillow` — runs on any Mac, Linux, or Windows box with a Python 3 environment and no HPC infrastructure.

## 📚 Citations

A publication describing bpp_identifier is in preparation. If you use the tool in your work in the meantime, please cite the underlying OpenCV framework along with any forthcoming Isom-Lab paper that features the tool.

1. Bradski G. *The OpenCV Library.* Dr. Dobb's Journal of Software Tools, 2000. [https://opencv.org](https://opencv.org)
2. Harris CR, Millman KJ, van der Walt SJ, et al. *Array programming with NumPy.* Nature. 2020 Sep;585(7825):357-362. doi: [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2). PMID: [32939066](https://pubmed.ncbi.nlm.nih.gov/32939066)
