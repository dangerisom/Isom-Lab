---
layout: default
title: Project vesicle_colocalization_quantifier
---

# Project vesicle_colocalization_quantifier

This project measures **two-channel pixel-level co-localization** of vesicles in microscopy images through an interactive, thresholded contour pipeline.

- [code](code/)
- [data](data/)

## 🔬 What vesicle_colocalization_quantifier Does

**vesicle_colocalization_quantifier** is a **Tkinter-based image-analysis application** that loads two microscopy channels (for example, a green vesicle marker and a red target marker), applies independent per-channel thresholds and contour detection (OpenCV), and computes the pixel-level **overlap** between the two. It reports the number of contours in each channel, the number of overlap regions, and the fraction of each channel that co-localizes with the other — then overlays the yellow intersection contours onto the composite image for visual confirmation. The GUI is designed for biologists who need reproducible vesicle co-localization counts without writing code.

## 🔍 Key Capabilities

- **Independent Per-Channel Thresholding**  
  Separate **low / mid / high** threshold sliders for **channel 1** (`t1`) and **channel 2** (`t2`), so a dim green channel and a bright red channel can each be dialed in without compromising the other.

- **Contour Area and Line-Width Bounds**  
  Configurable minimum and maximum contour areas and drawing line width prevent noise speckles from being counted as vesicles while still capturing small, legitimate puncta.

- **Pixel-Level Overlap Detection**  
  Builds boolean masks from each channel's contour set and computes their intersection, so "overlap" is measured at the pixel scale rather than at the coarse bounding-box level.

- **Yellow Overlap Contours**  
  Draws the intersection regions as **yellow** (BGR) contours on the composite overlay, giving immediate visual confirmation of which vesicles are scored as co-localized.

- **Quantitative Co-Localization Metrics**  
  Exports `total1`, `total2`, `overlap_count`, `frac1`, and `frac2` — the raw contour counts, the number of overlap regions, and the fraction of each channel that co-localizes with the other.

- **CSV Export**  
  Writes per-image metric rows to a CSV file, so a full directory of image pairs can be processed and aggregated without leaving the GUI.

## 🛠️ Additional Functionalities

- **Side-by-Side Display**  
  Shows the input channels next to the annotated overlap overlay so segmentation quality and overlap hits can be spot-checked before exporting.

- **Portable, Minimal Dependencies**  
  Requires only `cv2`, `numpy`, and the Python `csv` module — runs on any Mac, Linux, or Windows box with a Python 3 environment and no HPC infrastructure.

## 📚 Citations

A publication describing vesicle_colocalization_quantifier is in preparation. If you use the tool in your work in the meantime, please cite the underlying OpenCV framework along with any forthcoming Isom-Lab paper that features the tool.

1. Bradski G. *The OpenCV Library.* Dr. Dobb's Journal of Software Tools, 2000. [https://opencv.org](https://opencv.org)
2. Harris CR, Millman KJ, van der Walt SJ, et al. *Array programming with NumPy.* Nature. 2020 Sep;585(7825):357-362. doi: [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2). PMID: [32939066](https://pubmed.ncbi.nlm.nih.gov/32939066)
