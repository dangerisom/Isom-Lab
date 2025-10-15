#!/usr/bin/env python3
"""
Headless/cluster-friendly plotting pipeline (no CLI, no GUI).

Parses filenames ending with '-full.pdb' that look like:
  <rank>--<coverage>-<tm>--<anything>-<UNIPROT>-<anything>-full.pdb

Builds:
  dict_rank[uniprot]   = rank
  dict_scores[uniprot] = (coverage, tm)

Saves (in SAVE_DIR):
  results_scatter2d.png   # 2D scatter: coverage vs tm  (no labels)
  results_scatter3d.png   # 3D scatter: coverage, tm, rank  (no labels)
  results_contour2d.png   # 2D density contours ONLY (no overlay)
  results_contour3d.png   # 3D density surface ONLY  (no overlay)
  results_data.xlsx       # table of parsed rows

Optionally also saves per-rank plots under SAVE_DIR/by_rank/:
  rank_00012_scatter2d.png, rank_00012_contour2d.png, etc.

Notes:
- Uses SciPy KDE when dataset is small enough; otherwise a fast 2D histogram
  (with optional Gaussian blur) to keep runtime reasonable.
- Uses LogNorm with log-spaced contour levels for good contrast.
- Fixed x/y axis ranges for ALL plots via FIX_AXES/X/Y_MIN/MAX.
"""

# ======================== Configuration ========================
SOURCE_DIR = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.08.20.7tmp_sensitivity_test/2025.08.20.7tmps-rhod_mglu/2025.08.20.7tmps-rhod_mglu-analysis/2025.08.20.7tmps-rhod_mglu_results-filtered_hits/"   # e.g., "/path/to/input" or None for current working directory
SAVE_DIR   = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.08.20.7tmp_sensitivity_test/2025.08.20.7tmps-rhod_mglu/2025.08.20.7tmps-rhod_mglu-analysis/"   # e.g., "/path/to/output" or None to use SOURCE_DIR/cwd
OUTPUT_STEM = "rhod_mglu_results"

# Axis range controls (apply to ALL plots)
FIX_AXES = True
X_MIN, X_MAX = 0.45, 1.0   # default low/high for x
Y_MIN, Y_MAX = 0.45, 1.0   # default low/high for y

# Density grid / contour controls
GRID_SIZE        = 180        # 120–220 is typical. Bigger = smoother but slower.
KDE_BW           = None       # KDE bandwidth; None = Scott’s rule. Try 0.3–0.8 to sharpen.
KDE_MAX_POINTS   = 100_000    # Above this, auto-switch to fast histogram path
CHUNK_SIZE       = 50_000     # KDE evaluation chunk size for progress bar
SMOOTH_SIGMA_PX  = 1.5        # Gaussian blur sigma (grid cells) in fast path (0 to disable)
LEVELS_2D        = 20         # number of contour levels in 2D
LEVELS_3D        = 15         # number of contour levels projected on 3D base

# Per-rank plotting
PER_RANK_PLOTS       = False         # separate plots per rank
RANKS_TO_PLOT        = None         # e.g., [1,2,3] or None for all ranks present
MAX_RANKS_TO_PLOT    = None         # cap how many ranks to process; None = unlimited
PER_RANK_3D_CONTOUR  = True         # 3D density surface per rank
PER_RANK_3D_SCATTER  = False        # usually not useful (z=rank is constant)
RANK_SUBDIR_NAME     = "rhod_mglu_by_rank"
# ===============================================================

# Force a non-interactive backend for headless environments
import matplotlib
matplotlib.use("Agg")

from pathlib import Path
import re
import sys
from typing import Tuple, Dict, List, Optional
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, PowerNorm, ListedColormap
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 - needed for 3D projection
import pandas as pd
from tqdm import tqdm
from collections import defaultdict

# Optional SciPy imports
try:
    from scipy.stats import gaussian_kde
    SCIPY_KDE_AVAILABLE = True
except Exception:
    SCIPY_KDE_AVAILABLE = False

try:
    from scipy.ndimage import gaussian_filter
    SCIPY_ND_AVAILABLE = True
except Exception:
    SCIPY_ND_AVAILABLE = False


# ---------------------- Parsing ----------------------
def parse_filename(fname: str) -> Tuple[int, Tuple[float, float], str]:
    """
    Parse one filename using the specified rules.

    Returns
    -------
    rank : int
    scores : (coverage: float, tm: float)
    uniprot : str
    """
    if not fname.endswith("-full.pdb"):
        raise ValueError("Does not end with '-full.pdb'")

    stem = fname[:-9]  # strip '-full.pdb'
    parts = stem.split("--")
    if len(parts) < 3:
        raise ValueError("Expected at least 3 tokens when splitting on '--'")

    # 3a: rank
    rank_str = parts[0].strip()
    if not re.fullmatch(r"\d+", rank_str):
        raise ValueError(f"Rank is not an integer: {rank_str}")
    rank = int(rank_str)

    # 3b: coverage-tm
    score_part = parts[1].strip()
    score_bits = score_part.split("-")
    if len(score_bits) < 2:
        raise ValueError("Second token must contain 'coverage-tm'")
    coverage = float(score_bits[0])
    tm = float(score_bits[1])
    scores = (coverage, tm)

    # 3c / 3c.1: third token -> split on '-' and take SECOND token as UniProt
    token3 = parts[2].strip()
    t3_bits = token3.split("-")
    if len(t3_bits) < 2:
        raise ValueError("Third token must have at least two '-' parts to find UniProt ID")
    uniprot = t3_bits[1].strip()

    return rank, scores, uniprot


def collect_data(folder: Path):
    """
    Return:
      rows: List[Dict] with keys: filename, uniprot, rank, coverage, tm
      dict_rank:   {uniprot: rank} (lowest/best rank kept if duplicates)
      dict_scores: {uniprot: (coverage, tm)} aligned to best rank
    """
    rows: List[Dict] = []
    dict_rank: Dict[str, int] = {}
    dict_scores: Dict[str, Tuple[float, float]] = {}

    all_items = list(folder.iterdir())
    files = []
    for p in tqdm(all_items, desc="Scanning directory", unit="item", miniters=1, leave=False):
        if p.is_file() and p.name.endswith("-full.pdb"):
            files.append(p)

    if not files:
        print("No files ending with '-full.pdb' found.", file=sys.stderr)

    for p in tqdm(files, desc="Parsing files", unit="file", miniters=1):
        try:
            rank, (coverage, tm), uniprot = parse_filename(p.name)
        except Exception as e:
            tqdm.write(f"Skipping '{p.name}': {e}")
            continue

        row = {
            "filename": p.name,
            "uniprot": uniprot,
            "rank": rank,
            "coverage": coverage,
            "tm": tm,
        }
        rows.append(row)

        if (uniprot not in dict_rank) or (rank < dict_rank[uniprot]):
            dict_rank[uniprot] = rank
            dict_scores[uniprot] = (coverage, tm)

    return rows, dict_rank, dict_scores


# ---------------------- Axis helpers ----------------------
def _apply_axes_limits_2d():
    """Apply fixed axes to current 2D axes if FIX_AXES is True."""
    if FIX_AXES:
        plt.xlim(X_MIN, X_MAX)
        plt.ylim(Y_MIN, Y_MAX)

def _apply_axes_limits_3d(ax):
    """Apply fixed axes to given 3D axes if FIX_AXES is True."""
    if FIX_AXES:
        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(Y_MIN, Y_MAX)
        # z-axis (rank) left automatic


# ---------------------- Plots: scatter ----------------------
def make_scatter_2d(rows: List[Dict], out_png: Path):
    """2D scatter (coverage vs tm) without text labels."""
    if not rows:
        print("No rows to plot for 2D scatter.", file=sys.stderr)
        return
    x = [r["coverage"] for r in rows]
    y = [r["tm"] for r in rows]

    plt.figure(figsize=(7, 5))
    plt.scatter(x, y, s=40)
    _apply_axes_limits_2d()
    plt.title("Coverage vs TM")
    plt.xlabel("Coverage")
    plt.ylabel("TM score")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def make_scatter_3d(rows: List[Dict], out_png: Path):
    """3D scatter: x=coverage, y=tm, z=rank, without labels."""
    if not rows:
        print("No rows to plot for 3D scatter.", file=sys.stderr)
        return
    xs = [r["coverage"] for r in rows]
    ys = [r["tm"] for r in rows]
    zs = [r["rank"] for r in rows]

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(xs, ys, zs, s=30)
    _apply_axes_limits_3d(ax)
    ax.set_title("3D Scatter: Coverage (x), TM (y), Rank (z)")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("TM score")
    ax.set_zlabel("Rank")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


# ---------------------- Density grid helper ----------------------
def _density_grid_xy(x: np.ndarray, y: np.ndarray,
                     grid_size: int = GRID_SIZE,
                     kde_bw: Optional[float] = KDE_BW,
                     chunk_size: int = CHUNK_SIZE,
                     xlim: Optional[Tuple[float, float]] = None,
                     ylim: Optional[Tuple[float, float]] = None):
    """
    Build density grid Z over mesh (X, Y).
    - Uses SciPy gaussian_kde when available and dataset <= KDE_MAX_POINTS.
    - Otherwise uses a fast 2D histogram (optionally smoothed with Gaussian blur).
    Returns X, Y, Z with Z normalized to [0, 1].

    If xlim/ylim are provided, the grid exactly spans those ranges.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    n = x.size

    # Choose grid bounds
    if xlim is not None and ylim is not None:
        x_lo, x_hi = float(xlim[0]), float(xlim[1])
        y_lo, y_hi = float(ylim[0]), float(ylim[1])
    else:
        # derive from data with small padding
        if n == 0:
            x_lo, x_hi = 0.0, 1.0
            y_lo, y_hi = 0.0, 1.0
        else:
            dx = (x.max() - x.min())
            dy = (y.max() - y.min())
            pad_x = 0.05 * (dx if dx != 0 else 1.0)
            pad_y = 0.05 * (dy if dy != 0 else 1.0)
            x_lo, x_hi = x.min() - pad_x, x.max() + pad_x
            y_lo, y_hi = y.min() - pad_y, y.max() + pad_y

    # Build grid
    xi = np.linspace(x_lo, x_hi, grid_size)
    yi = np.linspace(y_lo, y_hi, grid_size)
    X, Y = np.meshgrid(xi, yi)

    # Empty input guard
    if n == 0:
        Z = np.zeros_like(X)
        return X, Y, Z

    use_kde = (SCIPY_KDE_AVAILABLE and n >= 2 and n <= KDE_MAX_POINTS)

    if use_kde:
        try:
            kde = gaussian_kde(np.vstack([x, y]), bw_method=kde_bw)
            pos = np.vstack([X.ravel(), Y.ravel()])
            Z = np.empty(pos.shape[1], dtype=float)
            for start in tqdm(range(0, pos.shape[1], chunk_size),
                              desc="Evaluating KDE grid", unit="pts", leave=False, miniters=1):
                end = min(start + chunk_size, pos.shape[1])
                Z[start:end] = kde(pos[:, start:end])
            Z = Z.reshape(X.shape)
        except Exception as e:
            tqdm.write(f"KDE fallback due to: {e}")
            use_kde = False

    if not use_kde:
        # Fast path: 2D histogram density
        H, xedges, yedges = np.histogram2d(
            x, y, bins=grid_size,
            range=[[x_lo, x_hi], [y_lo, y_hi]],
            density=True
        )
        # Center coords from bin edges
        Xc = 0.5 * (xedges[:-1] + xedges[1:])
        Yc = 0.5 * (yedges[:-1] + yedges[1:])
        X, Y = np.meshgrid(Xc, Yc)
        Z = H.T

        # Optional smoothing for nicer contours
        if SCIPY_ND_AVAILABLE and SMOOTH_SIGMA_PX and SMOOTH_SIGMA_PX > 0:
            Z = gaussian_filter(Z, sigma=SMOOTH_SIGMA_PX)

    # Normalize for consistent color scaling
    zmax = float(np.max(Z))
    if zmax > 0:
        Z = Z / zmax
    return X, Y, Z


# ---------------------- Plots: contours (separate) ----------------------
def make_contour_2d(rows: List[Dict], out_png: Path, grid_size: int = GRID_SIZE):
    """2D density contours ONLY (no scatter overlay), with lowest-density fill."""
    if not rows:
        print("No rows to plot for 2D contour.", file=sys.stderr)
        return
    x = np.array([r["coverage"] for r in rows])
    y = np.array([r["tm"] for r in rows])

    xlim = (X_MIN, X_MAX) if FIX_AXES else None
    ylim = (Y_MIN, Y_MAX) if FIX_AXES else None
    X, Y, Z = _density_grid_xy(x, y, grid_size=grid_size, xlim=xlim, ylim=ylim)

    # Positive range for LogNorm
    Zpos = Z[Z > 0]
    if Zpos.size == 0:
        print("All density values are zero; consider raising GRID_SIZE or SMOOTH_SIGMA_PX.", file=sys.stderr)
        return
    vmin = float(Zpos.min()); vmax = float(Zpos.max())
    eps = 1e-12
    vmin = max(vmin, eps)
    if vmax <= vmin:
        vmax = vmin * 10.0

    # Log-spaced levels (start at vmin)
    levels = np.geomspace(vmin, vmax, LEVELS_2D)

    # Key: values <= vmin go slightly BELOW vmin so 'under' is filled
    vfill = np.nextafter(vmin, 0.0)  # largest float < vmin
    Zplot = Z.copy()
    Zplot[Zplot <= vmin] = vfill

    # Build a cmap that supports 'under' across Matplotlib versions
    base_cmap = plt.get_cmap("viridis")
    try:
        # Matplotlib ≥ 3.3
        cmap = base_cmap.with_extremes(under=base_cmap(0.0))
    except AttributeError:
        # Older Matplotlib
        colors_array = base_cmap(np.linspace(0, 1, 256))
        cmap = ListedColormap(colors_array)
        cmap.set_under(base_cmap(0.0))

    fig, ax = plt.subplots(figsize=(7, 5))
    cntr = ax.contourf(
        X, Y, Zplot,
        levels=levels,
        cmap=cmap,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        extend="min"  # fill everything < vmin with 'under' color
    )
    ax.contour(
        X, Y, Zplot,
        levels=levels, colors="k", linewidths=0.3, alpha=0.35,
        norm=LogNorm(vmin=vmin, vmax=vmax)
    )

    if FIX_AXES:
        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_title("2D Density Contours: Coverage vs TM")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("TM score")
    ax.grid(True, linestyle="--", alpha=0.18)
    cbar = fig.colorbar(cntr); cbar.set_label("Relative point density")
    fig.tight_layout(); fig.savefig(out_png, dpi=200); plt.close(fig)


def make_contour_3d(rows: List[Dict], out_png: Path, grid_size: int = GRID_SIZE):
    """3D colored density surface ONLY (no scatter overlay), log-scaled.
       Zeros/negatives are clamped to vmin so base is lowest color."""
    if not rows:
        print("No rows to plot for 3D contour.", file=sys.stderr)
        return
    x = np.array([r["coverage"] for r in rows])
    y = np.array([r["tm"] for r in rows])

    xlim = (X_MIN, X_MAX) if FIX_AXES else None
    ylim = (Y_MIN, Y_MAX) if FIX_AXES else None
    X, Y, Z = _density_grid_xy(x, y, grid_size=grid_size, xlim=xlim, ylim=ylim)

    Zpos = Z[Z > 0]
    if Zpos.size == 0:
        print("All density values are zero; consider raising GRID_SIZE or SMOOTH_SIGMA_PX.", file=sys.stderr)
        return
    vmin = float(Zpos.min()); vmax = float(Zpos.max())
    eps = 1e-12
    vmin = max(vmin, eps)
    if vmax <= vmin:
        vmax = vmin * 10.0
    levels = np.geomspace(vmin, vmax, LEVELS_3D)

    # Clamp zeros/negatives up to vmin
    Zplot = Z.copy()
    Zplot[Zplot <= 0] = vmin

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(X, Y, Zplot, cmap="viridis", linewidth=0, antialiased=True,
                           norm=LogNorm(vmin=vmin, vmax=vmax), alpha=0.98)
    ax.contour(X, Y, Zplot, zdir="z", offset=vmin, levels=levels, cmap="viridis",
               norm=LogNorm(vmin=vmin, vmax=vmax), alpha=0.9)
    _apply_axes_limits_3d(ax)
    ax.set_title("3D Density Surface over Coverage–TM")
    ax.set_xlabel("Coverage"); ax.set_ylabel("TM score"); ax.set_zlabel("Relative density")
    fig.colorbar(surf, shrink=0.6, aspect=12, pad=0.08, label="Relative density")
    plt.tight_layout(); plt.savefig(out_png, dpi=200); plt.close()


# ---------------------- Excel export ----------------------
def save_excel(rows: List[Dict], out_xlsx: Path):
    if not rows:
        print("No rows to save to Excel.", file=sys.stderr)
        return
    df = pd.DataFrame(rows, columns=["uniprot", "rank", "coverage", "tm", "filename"])
    df = df.sort_values(by=["rank", "tm"], ascending=[True, False])
    df.to_excel(out_xlsx, index=False)


# ---------------------- Per-rank helpers ----------------------
def _group_rows_by_rank(rows: List[Dict]) -> Dict[int, List[Dict]]:
    buckets = defaultdict(list)
    for r in rows:
        buckets[r["rank"]].append(r)
    return dict(sorted(buckets.items(), key=lambda kv: kv[0]))


def make_per_rank_plots(rows: List[Dict], dst: Path):
    """Create separate scatter/contour plots for each rank (saved to by_rank/)."""
    rank_dir = dst / RANK_SUBDIR_NAME
    rank_dir.mkdir(parents=True, exist_ok=True)

    by_rank = _group_rows_by_rank(rows)
    all_ranks_sorted = list(by_rank.keys())

    # Which ranks?
    if RANKS_TO_PLOT is not None:
        target = [rk for rk in all_ranks_sorted if rk in set(RANKS_TO_PLOT)]
    else:
        target = all_ranks_sorted

    if MAX_RANKS_TO_PLOT is not None:
        target = target[:MAX_RANKS_TO_PLOT]

    for rk in tqdm(target, desc="Per-rank plots", unit="rank", miniters=1):
        grp = by_rank[rk]
        base = f"rank_{rk:05d}"
        f2d_scatter = rank_dir / f"{base}_scatter2d.png"
        f2d_contour = rank_dir / f"{base}_contour2d.png"
        f3d_scatter = rank_dir / f"{base}_scatter3d.png"
        f3d_contour = rank_dir / f"{base}_contour3d.png"

        try:
            make_scatter_2d(grp, f2d_scatter)
        except Exception as e:
            tqdm.write(f"[rank {rk}] scatter2d failed: {e}")

        try:
            make_contour_2d(grp, f2d_contour, grid_size=GRID_SIZE)
        except Exception as e:
            tqdm.write(f"[rank {rk}] contour2d failed: {e}")

        if PER_RANK_3D_SCATTER:
            try:
                make_scatter_3d(grp, f3d_scatter)
            except Exception as e:
                tqdm.write(f"[rank {rk}] scatter3d failed: {e}")

        if PER_RANK_3D_CONTOUR:
            try:
                make_contour_3d(grp, f3d_contour, grid_size=GRID_SIZE)
            except Exception as e:
                tqdm.write(f"[rank {rk}] contour3d failed: {e}")


# ---------------------- Main ----------------------
def main():
    # Validate axis config (avoid silent mistakes)
    if FIX_AXES:
        if not (X_MIN < X_MAX and Y_MIN < Y_MAX):
            raise ValueError(f"Invalid axis ranges: X_MIN({X_MIN})<X_MAX({X_MAX}), Y_MIN({Y_MIN})<Y_MAX({Y_MAX}) required.")

    src = Path(SOURCE_DIR).expanduser().resolve() if SOURCE_DIR else Path.cwd()
    dst = Path(SAVE_DIR).expanduser().resolve() if SAVE_DIR else src
    dst.mkdir(parents=True, exist_ok=True)

    print(f"Scanning (SOURCE_DIR): {src}")
    print(f"Saving   (SAVE_DIR)  : {dst}")
    if not SCIPY_KDE_AVAILABLE:
        print("Note: SciPy KDE not found — using fast histogram fallback.", file=sys.stderr)

    rows, dict_rank, dict_scores = collect_data(src)

    # # Log dictionaries (uncomment if useful)
    # print("\nDictionary (UniProt -> Rank):")
    # print(dict_rank)
    # print("\nDictionary (UniProt -> (Coverage, TM)):")
    # print(dict_scores)

    # Global outputs
    out_2d_scatter = dst / f"{OUTPUT_STEM}_scatter2d.png"
    out_3d_scatter = dst / f"{OUTPUT_STEM}_scatter3d.png"
    out_2d_contour = dst / f"{OUTPUT_STEM}_contour2d.png"
    out_3d_contour = dst / f"{OUTPUT_STEM}_contour3d.png"
    out_xlsx       = dst / f"{OUTPUT_STEM}_data.xlsx"

    make_scatter_2d(rows, out_2d_scatter)
    make_scatter_3d(rows, out_3d_scatter)
    make_contour_2d(rows, out_2d_contour, grid_size=GRID_SIZE)
    make_contour_3d(rows, out_3d_contour, grid_size=GRID_SIZE)
    save_excel(rows, out_xlsx)

    # Per-rank plots
    if PER_RANK_PLOTS:
        make_per_rank_plots(rows, dst)

    print("\nSaved:")
    print(f"  2D scatter : {out_2d_scatter}")
    print(f"  3D scatter : {out_3d_scatter}")
    print(f"  2D contour : {out_2d_contour}")
    print(f"  3D contour : {out_3d_contour}")
    print(f"  Excel data : {out_xlsx}")
    if PER_RANK_PLOTS:
        print(f"  Per-rank   : {dst / RANK_SUBDIR_NAME}/*")
    print("Done.")


if __name__ == "__main__":
    main()
