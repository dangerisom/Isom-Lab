---
layout: default
title: Project vesicle_triangulator
---
<!-- stager: preserve -->

# Project vesicle_triangulator

This project reconstructs **3D spatial relationships between vesicle transfer events** from 2D coordinate data using convex-hull / Delaunay triangulation, then writes the geometry out as both an annotated overlay image and an auto-formatted Excel workbook.

- [code](code/)
- [example_input](example_input/)
- [example_output](example_output/)

## 🔬 What vesicle_triangulator Does

**vesicle_triangulator** is a **Tkinter-based image-analysis application** that ingests CSV coordinate files plus their companion component images, builds a **3D convex hull and Delaunay triangulation** via the Isom-Lab `compGeometry` and `convexHull3D_2_1` modules, and emits the triangulated mesh as both a visualization overlay and a structured Excel report. It is designed for quantifying vesicle-transfer events between cells, where the raw data is a list of transfer-point coordinates that need to be placed into consistent 3D geometric context for downstream analysis.

## 🔍 Key Capabilities

- **CSV + Component-Image Ingestion**  
  Loads point coordinates from CSV files alongside their originating component images, so the triangulation is always anchored to the microscopy context in which the vesicle transfers were detected.

- **3D Convex Hull and Delaunay Triangulation**  
  Uses the Isom-Lab `compGeometry.Vertex` and `convexHull3D_2_1.convexHull3D` modules — the same production geometry stack used by pHinder — to compute the hull and Delaunay mesh over the vesicle transfer coordinates.

- **Annotated Overlay Image Output**  
  Renders the triangulation onto the component image so the geometric relationships between transfer events are visible directly on the microscopy substrate.

- **Auto-Fit Excel Workbook Export**  
  Writes the triangulation result — vertices, edges, and per-row attributes — into an `.xlsx` workbook whose column widths are auto-fitted via `openpyxl.utils.get_column_letter`, so the output is immediately readable without manual formatting.

- **Matplotlib-Backed Visualization**  
  Headless plotting of 2D/3D projections of the triangulation for figure generation and QA, using `matplotlib` alongside Pillow for raster compositing.

- **Pandas-Backed Data Handling**  
  Intermediate data are staged in `pandas` DataFrames so filters, sorts, and summary tables can be applied before the Excel export.

## 🛠️ Additional Functionalities

- **Shared Geometry Stack With pHinder**  
  Reuses the `compGeometry` and `convexHull3D_2_1` modules that power the pHinder side-chain-topology pipeline, so any improvements to the Isom-Lab hull engine propagate to vesicle triangulation automatically.

- **Single-File Tkinter GUI**  
  Bundled as a single-file `vesicle_transfer_triangulation_v10_1.py` application that runs with `python vesicle_transfer_triangulation_v10_1.py` on any Mac, Linux, or Windows box with a Python 3 environment.

## 📚 Citations

A publication describing vesicle_triangulator is in preparation. If you use the tool in your work in the meantime, please cite the underlying OpenCV, NumPy, and Matplotlib frameworks along with any forthcoming Isom-Lab paper that features the tool.

1. Bradski G. *The OpenCV Library.* Dr. Dobb's Journal of Software Tools, 2000. [https://opencv.org](https://opencv.org)
2. Harris CR, Millman KJ, van der Walt SJ, et al. *Array programming with NumPy.* Nature. 2020 Sep;585(7825):357-362. doi: [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2). PMID: [32939066](https://pubmed.ncbi.nlm.nih.gov/32939066)
3. Hunter JD. *Matplotlib: A 2D Graphics Environment.* Computing in Science & Engineering. 2007;9(3):90-95. doi: [10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55)
4. McKinney W. *Data Structures for Statistical Computing in Python.* Proceedings of the 9th Python in Science Conference, 2010;445:56-61. doi: [10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)
