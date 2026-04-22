# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.


#!/usr/bin/env python3
"""
UpSet diagrams + Excel (counts that match the UpSet plot), with simple tqdm bars.

Folder per domain:
parent/
  <domain>/
    *-analysis/
      *-filtered_hits/       (INCLUDE)
      *-filtered_hits_log/   (EXCLUDE)

File filter: only file names ending with "-full.pdb"
Filenames:   RANK--SCORE--TOKEN3[--...].ext

Install:
    pip install pandas openpyxl matplotlib upsetplot tqdm
"""

from pathlib import Path
from collections import defaultdict
from typing import Iterable, Dict, List, Tuple, Set

import os
import math
import warnings

import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships
from tqdm.auto import tqdm

# Quiet upsetplot/pandas warning under pandas 3.1
warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")

# ────────────────────────────
# Small helpers
# ────────────────────────────
def resolve_output_prefix(output_prefix: str | Path, parent: Path, default_stem: str = "upset_results") -> Path:
    p = Path(output_prefix) if output_prefix else Path("")
    if not str(p) or str(p) == ".":
        return parent / default_stem
    try:
        if p.exists() and p.is_dir():
            return p / default_stem
    except Exception:
        pass
    if str(p).endswith(("/", "\\")):
        return Path(str(p).rstrip("/\\")) / default_stem
    if p.suffix.lower() == ".xlsx":
        return p.with_suffix("")
    return p

def parse_filename_stem(stem: str):
    """
    Returns:
      (rank:int|None, score:str|None, token3:str|None, full_key:str|None)

    - token3   -> exactly the 3rd token after splitting on "--"
    - full_key -> token3 plus any remaining tokens (joined with "--")
    """
    parts = stem.split("--")
    if len(parts) < 3:
        return None, None, None, None
    rank_str, score_str = parts[0], parts[1]
    token3 = parts[2]
    full_key = "--".join(parts[2:])
    try:
        r = int(rank_str)
    except ValueError:
        return None, None, None, None
    return r, score_str, token3, full_key

# ────────────────────────────
# Folder walking (with "-full.pdb" filter)
# ────────────────────────────
SUFFIX_FILTER = "-full.pdb"  # change if needed

def iter_analysis_dirs(domain_path: Path):
    for p in domain_path.rglob("*-analysis"):
        if p.is_dir():
            yield p

def iter_filtered_hits_dirs(analysis_dir: Path):
    with os.scandir(analysis_dir) as it:
        for e in it:
            if e.is_dir(follow_symlinks=False):
                n = e.name
                if n.endswith("-filtered_hits_log"):
                    continue
                if n.endswith("-filtered_hits"):
                    yield Path(e.path)

def count_files_in_dir(dir_path: Path, name_endswith: str | tuple[str, ...] | None = None) -> int:
    with os.scandir(dir_path) as it:
        return sum(
            1 for d in it
            if d.is_file(follow_symlinks=False) and (name_endswith is None or d.name.endswith(name_endswith))
        )

def iter_files_in_dir(dir_path: Path, name_endswith: str | tuple[str, ...] | None = None):
    with os.scandir(dir_path) as it:
        for d in it:
            if d.is_file(follow_symlinks=False) and (name_endswith is None or d.name.endswith(name_endswith)):
                yield Path(d.path)

def scan_all_files(parent: Path) -> Tuple[List[Tuple[Path, str]], List[str]]:
    """Discover all files under <domain>/*-analysis/*-filtered_hits/* matching SUFFIX_FILTER."""
    domains = sorted([d for d in parent.iterdir() if d.is_dir()], key=lambda p: p.name)

    # total for a stable bar
    total = 0
    for domain in domains:
        for a in iter_analysis_dirs(domain):
            for h in iter_filtered_hits_dirs(a):
                total += count_files_in_dir(h, name_endswith=SUFFIX_FILTER)

    files_with_domains: List[Tuple[Path, str]] = []
    with tqdm(total=total, desc="Scanning", unit="file", leave=True) as bar:
        for domain in domains:
            dname = domain.name
            for a in iter_analysis_dirs(domain):
                for h in iter_filtered_hits_dirs(a):
                    for f in iter_files_in_dir(h, name_endswith=SUFFIX_FILTER):
                        files_with_domains.append((f, dname))
                        bar.update(1)
    return files_with_domains, [d.name for d in domains]

# ────────────────────────────
# Parsing + building sets (simple bar)
# ────────────────────────────
def parse_files(
    files_with_domains: List[Tuple[Path, str]],
    ranks: Iterable[int],
    postfix_every: int = 5000,   # update postfix every N files (0 to disable)
):
    ranks = tuple(int(r) for r in ranks)
    r2d2h: Dict[int, Dict[str, Set[str]]] = {r: defaultdict(set) for r in ranks}
    rows = []
    skipped = {"bad_pattern": 0, "bad_rank": 0, "rank_not_requested": 0}

    total = len(files_with_domains)
    with tqdm(total=total, desc="Parsing", unit="file", leave=True) as bar:
        for i, (f, domain) in enumerate(files_with_domains, start=1):
            r, score_str, token3, full_key = parse_filename_stem(f.stem)
            if r is None:
                if len(f.stem.split("--")) < 3:
                    skipped["bad_pattern"] += 1
                else:
                    skipped["bad_rank"] += 1
            elif r not in ranks:
                skipped["rank_not_requested"] += 1
            else:
                # IMPORTANT: use token3 for overlap logic
                r2d2h[r][domain].add(token3)
                rows.append({
                    "domain": domain,
                    "file": f.name,
                    "rank": r,
                    "score": score_str,
                    "token3": token3,
                    "hit_key_full": full_key,
                    "path": str(f),
                })

            bar.update(1)
            if postfix_every and (i % postfix_every == 0 or i == total):
                bar.set_postfix(skipped=sum(skipped.values()))

    return r2d2h, pd.DataFrame(rows), skipped

# ────────────────────────────
# Overlap + UpSet (aggregate duplicates; skip empties)
# ────────────────────────────
def compute_global_intersection_and_complements(domains: List[str], d2h: Dict[str, Set[str]]) -> Tuple[Set[str], Dict[str, Set[str]]]:
    sets = [d2h[d] for d in domains]
    if not sets:
        inter = set()
    elif len(sets) == 1:
        inter = set(sets[0])
    else:
        inter = set.intersection(*sets)
    comps = {d: d2h[d] - inter for d in domains}
    return inter, comps

def build_memberships_series(domains: List[str], d2h: Dict[str, Set[str]]) -> pd.Series:
    """Aggregated Series suitable for UpSet (no duplicate groups); empty Series if no hits."""
    if not domains:
        return pd.Series(dtype=int)
    all_hits = set().union(*[d2h.get(d, set()) for d in domains])
    if not all_hits:
        return pd.Series(dtype=int)

    memberships = []
    for hk in all_hits:
        present = tuple(d for d in domains if hk in d2h.get(d, set()))
        if present:
            memberships.append(present)
    if not memberships:
        return pd.Series(dtype=int)

    s = from_memberships(memberships)

    # Aggregate duplicate groups so UpSet can use subset_size="sum"
    if not s.index.is_unique:
        levels = list(range(s.index.nlevels))
        s = s.groupby(level=levels).sum()

    return s.astype(int, copy=False)

def build_upset_series_for_plot(domains: list[str],
                                d2h: dict[str, set[str]],
                                top_n: int,
                                min_size: int) -> pd.Series:
    """Return the exact Series the plot will use (apply min_size + top_n)."""
    s = build_memberships_series(domains, d2h)
    if s is None or getattr(s, "empty", True):
        return pd.Series(dtype=int)
    # filter by min_size (cardinality)
    if isinstance(min_size, int) and min_size > 1:
        s = s[s >= min_size]
    # sort and apply top_n
    s = s.sort_values(ascending=False)
    if isinstance(top_n, int) and top_n > 0 and len(s) > top_n:
        s = s.iloc[:top_n]
    return s

def upset_series_to_counts_df(s: pd.Series) -> pd.DataFrame:
    """Convert the UpSet Series to a tidy table: Combination | Count."""
    if s is None or getattr(s, "empty", True):
        return pd.DataFrame(columns=["Combination", "Count"])
    names = list(s.index.names)
    rows = []
    for key, val in s.items():
        # key is a tuple of booleans aligned with `names`
        if not isinstance(key, tuple):
            key = (key,)
        on = [name for name, present in zip(names, key) if bool(present)]
        combo = " & ".join(on) if on else "(none)"
        rows.append({"Combination": combo, "Count": int(val)})
    return pd.DataFrame(rows).sort_values("Count", ascending=False, kind="stable")

# ────────────────────────────
# Excel (summary + exact UpSet counts)
# ────────────────────────────
def save_excel_with_upset_counts(output_prefix: Path,
                                 ranks: Iterable[int],
                                 domains: list[str],
                                 r2d2h: dict[int, dict[str, set[str]]],
                                 top_n: int,
                                 min_size: int) -> Path:
    """
    Writes:
      - Summary (per-domain counts, global intersection size, complement size)
      - One sheet per rank: UpSet_Rank_<r> with EXACT plot combinations & counts
    """
    out_xlsx = Path(output_prefix).with_suffix(".xlsx")
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xlw:
        # Summary
        summary_rows = []
        for r in ranks:
            d2h = r2d2h[r]
            for d in domains:
                _ = d2h[d]  # ensure key
            inter, comps = compute_global_intersection_and_complements(domains, d2h)
            for d in domains:
                summary_rows.append({
                    "Rank": r,
                    "Domain": d,
                    "Domain_Set_Size": len(d2h[d]),
                    "Global_Intersection_Size": len(inter),
                    "Complement_Size": len(comps[d]),
                })
        pd.DataFrame(summary_rows).to_excel(xlw, sheet_name="Summary", index=False)

        # Exact plot counts per rank
        for r in ranks:
            s_plot = build_upset_series_for_plot(domains, r2d2h[r], top_n=top_n, min_size=min_size)
            df_counts = upset_series_to_counts_df(s_plot)
            sheet = f"UpSet_Rank_{r}"
            df_counts.to_excel(xlw, sheet_name=sheet, index=False)

    return out_xlsx


# ────────────────────────────
# CSV per rank (exact UpSet combinations & counts)
# ────────────────────────────
def save_rank_counts_csvs(output_prefix: Path,
                          ranks: Iterable[int],
                          domains: list[str],
                          r2d2h: dict[int, dict[str, set[str]]],
                          top_n: int,
                          min_size: int) -> list[Path]:
    """
    For each rank in `ranks`, write a CSV with the exact UpSet plot combinations and counts.
    Columns: Combination, Count
    Filename pattern: "<prefix>_rank<R>_upset_intersections.csv"
    """
    output_prefix = Path(output_prefix)
    written: list[Path] = []
    for r in ranks:
        d2h = r2d2h[r]
        # Build the exact series used for plotting (respects min_size and top_n)
        s_plot = build_upset_series_for_plot(domains, d2h, top_n=top_n, min_size=min_size)
        df_counts = upset_series_to_counts_df(s_plot)
        out_csv = output_prefix.with_name(f"{output_prefix.name}_rank{r}_upset_intersections.csv")
        df_counts.to_csv(out_csv, index=False)
        written.append(out_csv)
    return written

# ────────────────────────────
# Plotting
# ────────────────────────────
def save_upset_plot(output_prefix: Path, rank: int, domains: List[str],
                    d2h: Dict[str, Set[str]], top_n: int = 30, min_size: int = 1) -> Path | None:
    s_plot = build_upset_series_for_plot(domains, d2h, top_n=top_n, min_size=min_size)
    if s_plot.empty:
        return None
    upset = UpSet(
        s_plot,
        subset_size="sum",          # counts already aggregated
        sort_by="cardinality",
        sort_categories_by=None,
        intersection_plot_elements=10,
        min_subset_size=min_size,
    )
    plt.figure(figsize=(9, 6))
    upset.plot()
    plt.suptitle(f"UpSet Plot – Rank {rank}")
    out_png = Path(f"{output_prefix}_rank{rank}_upset.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return out_png

# ────────────────────────────
# Orchestration
# ────────────────────────────
def process_upset(
    parent_folder: str | Path,
    output_prefix: str | Path = ".",
    ranks: Iterable[int] = (1, 2, 3),
    top_n: int = 30,
    min_size: int = 1,
    write_raw_csv: bool = False,      # set True to write raw rows to chunked CSVs
    raw_chunk_rows: int = 1_000_000,
) -> dict:
    parent = Path(parent_folder).expanduser().resolve()
    safe_prefix = resolve_output_prefix(output_prefix, parent, default_stem="upset_results")

    files_with_domains, domain_names = scan_all_files(parent)
    r2d2h, raw_df, skipped = parse_files(files_with_domains, ranks)

    # ensure each domain key exists for each rank
    for r in ranks:
        for d in domain_names:
            _ = r2d2h[r][d]

    excel_path = save_excel_with_upset_counts(safe_prefix, ranks, domain_names, r2d2h,
                                              top_n=top_n, min_size=min_size)

    # optional raw CSV dump (kept off by default to avoid huge Excel)
    raw_csvs: List[Path] = []
    if write_raw_csv and not raw_df.empty:
        n = len(raw_df)
        chunks = math.ceil(n / raw_chunk_rows)
        for i in range(chunks):
            start, end = i * raw_chunk_rows, min((i + 1) * raw_chunk_rows, n)
            out_csv = Path(f"{safe_prefix}_raw_{i+1:02d}.csv")
            raw_df.iloc[start:end].to_csv(out_csv, index=False)
            raw_csvs.append(out_csv)

    images = {}
    for r in ranks:
        images[r] = save_upset_plot(safe_prefix, r, domain_names, r2d2h[r],
                                    top_n=top_n, min_size=min_size)
    # Write one CSV per rank with exact UpSet intersection counts
    csv_paths = save_rank_counts_csvs(safe_prefix, ranks, domain_names, r2d2h, top_n=top_n, min_size=min_size)

    return {"csv_paths": csv_paths,
            "domains": domain_names,
            "excel_path": excel_path,
            "images": images,
            "raw_csvs": raw_csvs,
            "skipped_counts": skipped,
            "num_files_scanned": len(files_with_domains),
            "num_rows_parsed": len(raw_df),
    }

# ────────────────────────────
# Example run (edit paths as needed)
# ────────────────────────────
if __name__ == "__main__":
    res = process_upset(
        parent_folder="/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.08.20.7tmp_sensitivity_test",
        output_prefix="/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.08.20.7tmp_sensitivity_test",
        ranks=(1, 2, 3),
        top_n=40,
        min_size=1,
        write_raw_csv=False,
    )
    print("Excel:", res["excel_path"])
    print("Images:", res["images"])
