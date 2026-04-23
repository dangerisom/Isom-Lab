# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UpSet plot for UNFILTERED hits — rock-solid version using from_indicators on the boolean matrix.

Outputs (SAVE_DIR or <PARENT_DIR>/aggregate_analysis):
  - aggregate_upset_membership_matrix.csv
  - aggregate_upset_plot.png / .svg         (UpSet grid)
  - aggregate_upset_top_intersections.png   (Top-K bar)
  - aggregate_upset_intersections.csv       (readable table of intersections)
  - aggregate_upset_per_domain_counts.png   (fallback when nothing to plot)
  - aggregate_upset_fallback_heatmap.png    (fallback visualization)
"""

# ===================== Configuration =====================
PARENT_DIR              = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.08.20.7tmp_sensitivity_test"

# Set to None to scan ALL domain folders under PARENT_DIR
# DOMAIN_DIR_PATTERN      = r"^2025\.08\.20\.7tmps\-rhod(?:_[^/]+)?$"
DOMAIN_DIR_PATTERN      = None


ANALYSIS_DIR_ENDSWITH   = "analysis"
CSV_FILENAME_SUFFIX     = "_results.csv"

SAVE_DIR                = None

# Plot tuning
UPSET_FIGSIZE           = (12, 7)
TOP_INTERSECTIONS_K     = 20

# UniProt parsing options
STRIP_ISOFORM_SUFFIX    = True
WRITE_DEBUG_LISTS       = False
# =========================================================

import os, csv, re, warnings
from pathlib import Path
from typing import List, Tuple, Dict, Set, Optional
from collections import Counter

import numpy as np
import pandas as pd

from tqdm import tqdm

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

# Quiet a noisy futurewarning from upsetplot+pandas inside their internals.
try:
    warnings.filterwarnings(
        "ignore",
        message=r"Downcasting object dtype arrays on \.fillna",
        module=r"upsetplot\.data"
    )
    pd.set_option('future.no_silent_downcasting', True)
except Exception:
    pass

# upsetplot imports
try:
    from upsetplot import UpSet, from_indicators
    _HAS_FROM_INDICATORS = True
except Exception:
    _HAS_FROM_INDICATORS = False
    try:
        from upsetplot import UpSet
    except Exception:
        UpSet = None  # type: ignore


# ---------- Utilities ----------
def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _filter_domain_dirs(all_dirs: List[Path]) -> List[Path]:
    if DOMAIN_DIR_PATTERN in (None, "", False):
        return all_dirs
    pat = re.compile(DOMAIN_DIR_PATTERN, re.IGNORECASE)
    return [d for d in all_dirs if pat.search(d.name) is not None]


def _list_domain_csvs(parent: Path) -> List[Tuple[str, Path]]:
    pairs: List[Tuple[str, Path]] = []

    all_domain_dirs: List[Path] = sorted([d for d in parent.iterdir() if d.is_dir()])
    if not all_domain_dirs:
        raise FileNotFoundError(f"No domain subdirectories found in {parent}")

    domain_dirs: List[Path] = _filter_domain_dirs(all_domain_dirs)
    if not domain_dirs:
        raise FileNotFoundError(
            f"No domain folders matched DOMAIN_DIR_PATTERN={DOMAIN_DIR_PATTERN!r} under {parent}"
        )

    print(f"Discovered {len(all_domain_dirs)} domain folders; "
          f"{len(domain_dirs)} matched filter {DOMAIN_DIR_PATTERN!r}.")
    print("Matched domain folders:")
    for d in domain_dirs:
        print(f"  - {d.name}")

    totals = Counter()
    t_domains = tqdm(domain_dirs, desc="Scanning domains", unit="domain", leave=True, dynamic_ncols=True)
    ends = ANALYSIS_DIR_ENDSWITH.lower()
    csv_suffix = CSV_FILENAME_SUFFIX

    for domain_dir in t_domains:
        domain = domain_dir.name
        if domain == "aggregate_analysis":
            continue
        t_domains.set_postfix_str(domain, refresh=False)

        # 1) Immediate children check for *analysis_test
        analysis_dirs_immediate: List[Path] = []
        try:
            for child in domain_dir.iterdir():
                if child.is_dir() and child.name.lower().endswith(ends):
                    analysis_dirs_immediate.append(child)
        except PermissionError:
            analysis_dirs_immediate = []

        csvs_found: List[Path] = []
        for adir in analysis_dirs_immediate:
            for f in adir.glob(f"*{csv_suffix}"):
                if f.is_file():
                    csvs_found.append(f)

        used_fallback = False
        analysis_dirs_count = len(analysis_dirs_immediate)

        # 2) If none, single streaming os.walk
        if not csvs_found:
            preferred_hits: List[Path] = []
            fallback_hits: List[Path] = []
            for root, dirs, files in os.walk(domain_dir, topdown=True):
                base = os.path.basename(root).lower()
                in_analysis = base.endswith(ends)

                if in_analysis:
                    for name in files:
                        if name.endswith(csv_suffix):
                            preferred_hits.append(Path(root) / name)
                for name in files:
                    if name.endswith(csv_suffix):
                        fallback_hits.append(Path(root) / name)

            if preferred_hits:
                csvs_found = sorted(set(preferred_hits))
                analysis_dirs_count = len({p.parent for p in preferred_hits})
            else:
                used_fallback = True
                csvs_found = sorted(set(fallback_hits))

        totals["domains"] += 1
        totals["analysis_dirs"] += analysis_dirs_count
        totals["csvs"] += len(csvs_found)
        if used_fallback and csvs_found:
            totals["fallback_domains_with_csvs"] += 1
        if used_fallback and not csvs_found:
            totals["fallback_domains_without_csvs"] += 1

        tqdm.write(f"[{domain}] analysis_dirs={analysis_dirs_count} "
                   f"csvs_found={len(csvs_found)} "
                   f"{'(fallback)' if used_fallback else ''}")

        for f in csvs_found:
            pairs.append((domain, f))

    tqdm.write("—" * 72)
    tqdm.write(f"Discovery summary: domains={totals['domains']}, "
               f"analysis_dirs={totals['analysis_dirs']}, csvs={totals['csvs']}, "
               f"fallback_with_csvs={totals['fallback_domains_with_csvs']}, "
               f"fallback_without_csvs={totals['fallback_domains_without_csvs']}")
    tqdm.write("—" * 72)

    return pairs


# ---------- UniProt parsing (robust) ----------
_RE_AF   = re.compile(r"AF-([A-Za-z0-9]{6,10})-", re.IGNORECASE)
_RE_ACC6 = re.compile(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]{5})", re.IGNORECASE)
_RE_ACC10= re.compile(r"[A-NR-Z][0-9][A-Z0-9]{8}", re.IGNORECASE)

def _normalize_uniprot(uid: str) -> str:
    u = uid.strip().upper()
    if STRIP_ISOFORM_SUFFIX:
        m = re.match(r"^([A-Z0-9]{6,10})(?:-\d+)?$", u)
        if m:
            return m.group(1)
    return u

def _parse_uniprot_from_row_fields(fields: List[str]) -> Optional[str]:
    if len(fields) >= 4:
        f = fields[3].strip()
        bits = f.split("-")
        if len(bits) >= 2:
            return _normalize_uniprot(bits[1])

    for f in fields:
        m = _RE_AF.search(f)
        if m:
            return _normalize_uniprot(m.group(1))

    for f in fields:
        m = _RE_ACC10.search(f)
        if m:
            return _normalize_uniprot(m.group(0))

    for f in fields:
        m = _RE_ACC6.search(f)
        if m:
            return _normalize_uniprot(m.group(0))

    return None


# ---------- Core ----------
def collect_domain_sets(parent: Path) -> Dict[str, Set[str]]:
    domain_sets: Dict[str, Set[str]] = {}
    csvs = _list_domain_csvs(parent)

    if not csvs:
        raise FileNotFoundError(
            f"No CSVs ending with '{CSV_FILENAME_SUFFIX}' under {parent}/<domain>/..."
        )

    for domain, csv_path in tqdm(csvs, desc="Reading CSVs", unit="file", leave=True, dynamic_ncols=True):
        ids = domain_sets.setdefault(domain, set())
        with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as fh:
            reader = csv.reader(fh)
            first = True
            for row in reader:
                if not row:
                    continue
                if first:
                    first = False
                    hdr = ",".join(row)
                    if any(h in hdr for h in ("TM score", "Query", "PDB", "Target", "filename", "header")):
                        continue  # header
                uid = _parse_uniprot_from_row_fields(row)
                if uid:
                    ids.add(uid)

    if WRITE_DEBUG_LISTS:
        dbg_dir = (Path(SAVE_DIR).expanduser().resolve() if SAVE_DIR else Path(PARENT_DIR) / "aggregate_analysis")
        _ensure_dir(dbg_dir)
        for d, s in domain_sets.items():
            with (dbg_dir / f"debug_ids_{d}.txt").open("w") as f:
                for uid in sorted(s):
                    f.write(uid + "\n")

    domain_sets = {d: s for d, s in domain_sets.items() if s}
    if not domain_sets:
        raise RuntimeError("No UniProt IDs parsed. Check CSV format or parsing rule.")
    return domain_sets


def build_membership_matrix(domain_sets: Dict[str, Set[str]]) -> pd.DataFrame:
    all_ids = sorted(set().union(*domain_sets.values()))
    domains = sorted(domain_sets.keys())
    df = pd.DataFrame(0, index=all_ids, columns=domains, dtype=int)
    for d in domains:
        if domain_sets[d]:
            df.loc[list(domain_sets[d]), d] = 1
    return df


# ---------- Plot helpers ----------
def _labels_for_intersections(series_counts: pd.Series, max_len: int = 80) -> list:
    labels = []
    names = list(series_counts.index.names)
    if any(n is None for n in names):
        names = [f"Set{i+1}" for i in range(series_counts.index.nlevels)]
    for tup in series_counts.index:
        active = [name for name, is_on in zip(names, tup) if bool(is_on)]
        lbl = " + ".join(active) if active else "(none)"
        labels.append(lbl if len(lbl) <= max_len else lbl[:max_len-1] + "…")
    return labels


def _to_numeric_series_1d(obj) -> pd.Series:
    if isinstance(obj, pd.DataFrame):
        if obj.shape[1] == 1:
            s = obj.iloc[:, 0]
        elif "count" in obj.columns:
            s = obj["count"]
        else:
            num = obj.select_dtypes(include=[np.number])
            s = num.sum(axis=1) if not num.empty else pd.to_numeric(obj.iloc[:, 0], errors="coerce")
    else:
        s = obj
    if not isinstance(s, pd.Series):
        s = pd.Series(s)
    return pd.to_numeric(s, errors="coerce").fillna(0).astype("int64")


def _plot_top_intersections_bar_from_series(series_counts, out_png: Path, top_k: int = TOP_INTERSECTIONS_K) -> Optional[Path]:
    series_counts = _to_numeric_series_1d(series_counts)
    nz = series_counts[series_counts.gt(0)].sort_values(ascending=False)
    if nz.empty:
        return None
    top = nz.head(top_k)
    labels = _labels_for_intersections(top)
    width = min(24, max(8, 0.4 * len(top)))
    fig, ax = plt.subplots(figsize=(width, 6))
    ax.bar(np.arange(len(top)), top.values)
    ax.set_xticks(np.arange(len(top)))
    ax.set_xticklabels(labels, rotation=90, ha="center")
    ax.set_ylabel("Intersection size (IDs)")
    ax.set_title(f"Top {min(top_k, len(nz))} intersections")
    fig.tight_layout()
    bar_png = out_png.with_name(out_png.stem.replace("upset", "top_intersections") + out_png.suffix)
    fig.savefig(bar_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return bar_png


def plot_fallback_heatmap(df_bool: pd.DataFrame, out_png: Path) -> None:
    if df_bool.empty:
        return
    MAX_IDS = 1000
    df_small = df_bool.iloc[:MAX_IDS]
    fig, ax = plt.subplots(figsize=(max(8, len(df_small.columns) * 0.5), 8))
    im = ax.imshow(df_small.values, aspect="auto", interpolation="nearest", cmap="Greys")
    ax.set_yticks(np.arange(len(df_small.index)))
    if len(df_small.index) > 40:
        ax.set_yticks(np.linspace(0, len(df_small.index) - 1, 40).astype(int))
        ax.set_yticklabels(df_small.index[ax.get_yticks()])
    else:
        ax.set_yticklabels(df_small.index)
    ax.set_xticks(np.arange(len(df_small.columns)))
    ax.set_xticklabels(df_small.columns, rotation=90)
    ax.set_title("Membership heatmap (fallback)")
    ax.set_xlabel("Domains")
    ax.set_ylabel("UniProt IDs")
    fig.colorbar(im, ax=ax, label="Membership (1 = present)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ---------- UpSet pipeline (robust) ----------
def build_upset_series_from_bool_df(bool_df: pd.DataFrame) -> pd.Series:
    """
    Version-proof: derive intersection counts by grouping across all boolean columns.
    Returns a 1-D int64 Series with a boolean MultiIndex (one level per domain).
    """
    # Ensure strict boolean indicator matrix
    dfb = bool_df.astype(bool)

    # Group by all domain columns => counts per True/False combination
    cols = list(dfb.columns)
    s = dfb.groupby(cols, dropna=False).size()
    s.index = s.index.set_names(cols)

    # Drop the 'all-False' bucket (items that belong to no domain), ignore if absent
    all_false = tuple([False] * len(cols))
    s = s.drop(index=all_false, errors="ignore")

    # Keep only nonzero intersections, make sure dtype is numeric and sorted
    s = s[s > 0].astype("int64").sort_values(ascending=False)
    return s


def plot_upset_from_bool_df(domain_sets: Dict[str, Set[str]], out_png: Path, df_bool: pd.DataFrame) -> None:
    bool_df = df_bool.astype(bool).copy()
    if bool_df.empty or bool_df.shape[1] == 0:
        print("[warn] Empty membership matrix; writing heatmap fallback only.")
        plot_fallback_heatmap(bool_df, out_png.with_name(out_png.stem.replace("upset", "fallback_heatmap") + out_png.suffix))
        return

    # Build counts series from the boolean indicator matrix
    series = build_upset_series_from_bool_df(bool_df)
    series = series.sort_values(ascending=False)

    # Save readable intersections CSV
    names = list(series.index.names)
    rows = []
    for key, val in series.items():
        if not isinstance(key, tuple):
            key = (key,)
        on = [name for name, present in zip(names, key) if bool(present)]
        combo = " & ".join(on) if on else "(none)"
        rows.append({"Combination": combo, "Count": int(val)})
    inter_csv = out_png.with_name(out_png.stem + "_intersections.csv")
    pd.DataFrame(rows).sort_values("Count", ascending=False, kind="stable").to_csv(inter_csv, index=False)
    print(f"[diag] Intersections table: {inter_csv} (combos={len(series)}, total={int(series.sum())})")

    # Top-K bar (always useful)
    _plot_top_intersections_bar_from_series(series, out_png, top_k=TOP_INTERSECTIONS_K)

    if series.empty:
        print("[warn] No nonzero intersections; writing fallbacks.")
        dom_sizes = pd.Series({d: len(s) for d, s in domain_sets.items()}).sort_values(ascending=False)
        if not dom_sizes.empty:
            fig, ax = plt.subplots(figsize=(max(8, len(dom_sizes) * 0.6), 5))
            ax.bar(dom_sizes.index, dom_sizes.values)
            ax.set_ylabel("Count of UniProt IDs"); ax.set_xlabel("Domains")
            ax.set_title("Per-domain counts (no overlaps to plot)")
            plt.xticks(rotation=90); plt.tight_layout()
            fig.savefig(out_png.with_name(out_png.stem.replace("upset", "per_domain_counts") + out_png.suffix),
                        dpi=200, bbox_inches="tight")
            plt.close(fig)
        plot_fallback_heatmap(bool_df, out_png.with_name(out_png.stem.replace("upset", "fallback_heatmap") + out_png.suffix))
        return

    # Cap to Top-K for readability in the grid as well
    series_vis = series.head(TOP_INTERSECTIONS_K)

    if UpSet is None:
        # ... (unchanged)
        return

    # Actual UpSet grid
    try:
        # upset = UpSet(
        #     series_vis,
        #     sort_by="cardinality",
        #     show_counts="%d",
        #     subset_size="sum",      # <-- tell UpSet to use our aggregated counts
        # )
        upset = UpSet(
            series_vis,
            sort_by="cardinality",
            show_counts=False,   # hide numeric labels
            subset_size="sum",
        )
        upset.plot()
        fig = plt.gcf()
        fig.set_size_inches(*UPSET_FIGSIZE)
        plt.suptitle("UpSet: UniProt ID overlaps across query domains", y=0.98, fontsize=14)
        plt.tight_layout()
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        try:
            fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight")
        except Exception:
            pass
        plt.close(fig)
        print(f"Saved UpSet grid: {out_png}")
    except Exception as e:
        print(f"[warn] UpSet failed: {e}; writing fallbacks.")
        dom_sizes = pd.Series({d: len(s) for d, s in domain_sets.items()}).sort_values(ascending=False)
        if not dom_sizes.empty:
            fig, ax = plt.subplots(figsize=(max(8, len(dom_sizes) * 0.6), 5))
            ax.bar(dom_sizes.index, dom_sizes.values)
            ax.set_ylabel("Count of UniProt IDs"); ax.set_xlabel("Domains")
            ax.set_title("Per-domain counts (plot fallback)")
            plt.xticks(rotation=90); plt.tight_layout()
            fig.savefig(out_png.with_name(out_png.stem.replace("upset", "per_domain_counts") + out_png.suffix),
                        dpi=200, bbox_inches="tight")
            plt.close(fig)
        plot_fallback_heatmap(bool_df, out_png.with_name(out_png.stem.replace("upset", "fallback_heatmap") + out_png.suffix))


# ---------------------------------- Main ------------------------------------- #
def main() -> None:
    parent = Path(PARENT_DIR).expanduser().resolve()
    out_dir = Path(SAVE_DIR).expanduser().resolve() if SAVE_DIR else (parent / "aggregate_analysis")
    _ensure_dir(out_dir)

    print(f"Parent directory      : {parent}")
    print(f"Output directory      : {out_dir}")
    print(f"Domain filter (regex) : {DOMAIN_DIR_PATTERN!r}")
    print(f"Analysis dir endswith : '{ANALYSIS_DIR_ENDSWITH}'")
    print(f"CSV filename suffix   : '{CSV_FILENAME_SUFFIX}'")

    # 1) Gather sets per domain
    domain_sets = collect_domain_sets(parent)
    print(f"Domains with ≥1 ID: {len(domain_sets)}")
    for d, s in sorted(domain_sets.items()):
        print(f"  - {d}: {len(s)} unique UniProt IDs")

    # 2) Build membership matrix (IDs x Domains) and save
    df = build_membership_matrix(domain_sets)
    matrix_csv = out_dir / "aggregate_upset_membership_matrix.csv"
    df.to_csv(matrix_csv, index=True)
    print(f"Saved membership matrix: {matrix_csv}")

    # 3) UpSet plot built from the boolean matrix (no API quirks)
    upset_png = out_dir / "aggregate_upset_plot.png"
    plot_upset_from_bool_df(domain_sets, upset_png, df.astype(bool))
    print(f"Attempted UpSet plot   : {upset_png}")
    print("Done.")


if __name__ == "__main__":
    main()
