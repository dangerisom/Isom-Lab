# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.


#!/usr/bin/env python3
"""
Go-forward network builder from thresholded BLAST similarity graph.

Inputs (from your previous network builder):
  network_out/nodes.tsv.gz
      (uniprot_id, superkingdom, gene, organism, rank, degree_unique)
  network_out/edges_unique.tsv.gz
      Supported formats:
        • Legacy 3 columns:  src_uid, dst_uid, pident               (evalue/bitscore treated as missing)
        • Extended 6 columns: src_uid, dst_uid, pident, evalue, length, bitscore

What this script does
---------------------
1) Rebuild a weighted, undirected graph from edges_unique.tsv.gz,
   keeping only edges with:
     • max_pident >= IDENTITY_THRESHOLD, and
     • a significance rule based on FILTER_MODE:
         - FILTER_MODE == "evalue": require finite evalue <= EVALUE_MAX
         - FILTER_MODE == "bitscore": require finite bitscore >= BITSCORE_MIN
   Both weight (pident), evalue, and bitscore are carried through traversal and outputs.

2) Choose a start node:
     - Prefer the highest-degree node whose superkingdom matches START_SUPERKINGDOM
       (after normalization), if any such node exists in the thresholded graph.
     - Otherwise, fall back to the highest-degree node whose superkingdom is one of
       ALLOWED_SUPERKINGDOMS (Archaea/Bacteria/Eukaryota).

3) Build a "go-forward" network by BFS from the start node, obeying:
     - Once a node is encountered, it cannot be re-used (no revisits => no rings).
     - Only traverse edges that meet the %identity threshold and the selected significance rule.
     - Neighbor expansion order per node is:
         (highest pident) → (lowest evalue; NaN/absent treated worst) → (lexicographic tiebreaker)

4) Write outputs that preserve connections (folder & filenames include thresholds):
     go_forward_out_id<ID>_<KVAL>[_<RUN_SUFFIX>]/go_forward_nodes_id<ID>_<KVAL>[_<RUN_SUFFIX>].tsv.gz
         (node, parent, weight_to_parent, depth, component_id, attrs...)
     go_forward_out_id<ID>_<KVAL>[_<RUN_SUFFIX>]/go_forward_edges_id<ID>_<KVAL>[_<RUN_SUFFIX>].tsv.gz
         (parent, child, weight, evalue, bitscore, component_id)
     go_forward_out_id<ID>_<KVAL>[_<RUN_SUFFIX>]/go_forward_summary_id<ID>_<KVAL>[_<RUN_SUFFIX>].txt

     Where:
       • <ID>   = IDENTITY_THRESHOLD as an integer (e.g., 30 for 30.0)
       • <KVAL> = "e<EVALUE_MAX>" if FILTER_MODE == "evalue" (e.g., e0.1),
                  or "b<BITSCORE_MIN>" if FILTER_MODE == "bitscore" (e.g., b60)
       • RUN_SUFFIX is an optional user-provided tag to uniquely identify the run

Stdlib only. Large-file friendly (streams inputs).
"""

from __future__ import annotations
import gzip
from dataclasses import dataclass
from pathlib import Path
from collections import deque, defaultdict
from typing import Dict, List, Tuple, Iterable, Optional
import math

def _format_run_suffix(
    id_thr: float,
    emax: float,
    run_suffix: str,
    filter_mode: str = "evalue",
    bsmin: float = float("nan"),
) -> str:
    # choose tag based on filter mode
    if (filter_mode or "").lower() == "bitscore":
        k_str = (f"b{bsmin:g}" if math.isfinite(bsmin) else "bNaN")
    else:
        k_str = ("einf" if math.isinf(emax) else f"e{emax:g}")
    safe = "".join(ch for ch in run_suffix if ch.isalnum() or ch in ("-", "_", "."))
    base = f"id{int(id_thr)}_{k_str}"
    return base if not safe else base + f"_{safe}"


# =========================== CONFIG ===========================
# Path to the folder that contains network_out/* from the previous step
# NETWORK_DIR = Path("/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/2023.11.12.7tmps-analysis_final/2023.11.12.trim_blastp_sequences_final_nr_id0p50_c0p50/blastp_out/network_out_fast")
# NETWORK_DIR = Path("/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/2023.11.12.7tmps-analysis_final/ABE_reps/blastp_out/network_out_fast")
NETWORK_DIR = Path("/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/2023.11.12.7tmps-analysis_final/trim_blastp_sequences_human_gpcrs/blastp_out/network_out_fast")


# Identity threshold (percent identity) for including edges in the go-forward network
IDENTITY_THRESHOLD: float = 30.0
# Maximum allowed e-value for edges during traversal (set to inf to disable)
# EVALUE_MAX: float = float('inf')
EVALUE_MAX: float = float(0.1)
# Optional user-provided suffix to uniquely identify this run's folder/files (safe chars only)
# RUN_SUFFIX: str = "abe_cu"
RUN_SUFFIX: str = "hsapiens_gpcrs_cu"

# Filtering mode for traversal: "evalue" or "bitscore"
FILTER_MODE: str = "bitscore"
# Minimum bitscore for edges during traversal (used when FILTER_MODE == "bitscore")
BITSCORE_MIN: float = 30.0

# Start from this superkingdom if possible (preferred seed choice).
# Use one of: "archaea", "bacteria", "eukaryota" (case-insensitive), or None to disable preference.
# START_SUPERKINGDOM: Optional[str] = "archaea"
START_SUPERKINGDOM: Optional[str] = "bacteria"

# Which superkingdom labels qualify as "labeled starts" for fallback selection.
# These are normalized via _norm_sk before comparison.
ALLOWED_SUPERKINGDOMS = {"archaea", "bacteria", "eukaryota"}

# Where to write the outputs
OUT_DIR = NETWORK_DIR / (f"go_forward_out_" + _format_run_suffix(IDENTITY_THRESHOLD, EVALUE_MAX, RUN_SUFFIX, FILTER_MODE, BITSCORE_MIN))
OVERWRITE_OUTPUT = True  # if True and OUT_DIR exists, it will be removed and recreated
# =============================================================


@dataclass
class NodeAttr:
    superkingdom: str
    gene: str
    organism: str
    rank: str
    degree_unique: int


def _ensure_outdir(out: Path, overwrite: bool = False) -> None:
    if out.exists() and overwrite:
        # careful, delete only the leaf directory
        import shutil
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)


def _read_nodes(nodes_path: Path) -> Dict[str, NodeAttr]:
    """
    nodes.tsv.gz format (tab):
      uniprot_id  superkingdom  gene  organism  rank  degree_unique
    """
    nodes: Dict[str, NodeAttr] = {}
    with gzip.open(nodes_path, "rt", encoding="utf-8") as fh:
        header = next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            uid, sk, gene, org, rank, deg = parts[:6]
            try:
                d = int(deg)
            except Exception:
                d = 0
            nodes[uid] = NodeAttr(superkingdom=sk, gene=gene, organism=org, rank=rank, degree_unique=d)
    return nodes


def _read_weighted_edges(edges_path: Path, min_identity: float, max_evalue: float = float('inf')) -> Dict[str, List[Tuple[str, float, float, float]]]:
    """
    edges_unique.tsv.gz expected formats (tab):
      legacy: src_uid  dst_uid  max_pident
      new   : src_uid  dst_uid  pident  evalue  length  bitscore
    Returns adjacency dict filtered by min_identity (undirected), as (neighbor, pident, evalue).
    If evalue is not present, it is set to float("nan").
    """
    adj: Dict[str, List[Tuple[str, float, float]]] = defaultdict(list)
    with gzip.open(edges_path, "rt", encoding="utf-8") as fh:
        header = next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            a, b = parts[0], parts[1]
            try:
                if len(parts) >= 6:
                    w = float(parts[2])
                    ev = float("inf") if parts[3] == "inf" else float(parts[3])
                    bs = float(parts[5])
                else:
                    w = float(parts[2])
                    ev = float("nan")
                    bs = float("nan")
            except Exception:
                continue
            if w >= min_identity and (
                (FILTER_MODE.lower() == 'evalue' and (math.isinf(max_evalue) or (math.isfinite(ev) and ev <= max_evalue))) or
                (FILTER_MODE.lower() == 'bitscore' and (math.isfinite(bs) and bs >= BITSCORE_MIN))
            ):
                adj[a].append((b, w, ev, bs))
                adj[b].append((a, w, ev, bs))
    return adj



def _norm_sk(sk: str) -> str:
    s = (sk or "").strip().lower()
    # harmonize common variants
    if s.startswith("arch"):                 # archaea / archaeal
        return "archaea"
    if s.startswith("bact"):                 # bacteria / bacterial
        return "bacteria"
    if s.startswith(("euk", "euka", "eukar", "eukary", "eukaryot")):  # eukaryota / eukaryotic / eukaryote
        return "eukaryota"
    return s


def _choose_start_node(adj: Dict[str, List[Tuple[str, float, float]]],
                       attrs: Dict[str, NodeAttr]) -> Optional[str]:
    """
    Pick highest-degree node (in the *thresholded* graph).
    Preference order:
      1) If START_SUPERKINGDOM is set, pick the highest-degree node normalized to that SK.
      2) Otherwise (or if none found), pick the highest-degree node among ALLOWED_SUPERKINGDOMS.
    """
    def deg(u: str) -> int:
        return len(adj.get(u, ()))
    def sk_of(u: str) -> str:
        a = attrs.get(u)
        return _norm_sk(a.superkingdom if a else "")

    # 1) Prefer a specific superkingdom if requested
    if START_SUPERKINGDOM:
        target = _norm_sk(START_SUPERKINGDOM)
        candidates = [u for u in adj.keys() if sk_of(u) == target]
        if candidates:
            return max(candidates, key=lambda u: (deg(u), u))

    # 2) Fallback to any allowed labeled SK
    candidates = [u for u in adj.keys() if sk_of(u) in ALLOWED_SUPERKINGDOMS]
    return max(candidates, key=lambda u: (deg(u), u)) if candidates else None


def _bfs_go_forward(adj: Dict[str, List[Tuple[str, float, float]]],
                    root: str) -> Tuple[List[Tuple[str, Optional[str], float, int, int]],
                                        List[Tuple[str, str, float, float, int]]]:
    """
    BFS from root on an undirected, weighted graph to produce a directed, acyclic
    "go forward" tree (no rings). Returns:

      nodes_out: [(node, parent, weight_to_parent, depth, component_id)]
      edges_out: [(parent, child, weight, component_id)]

    component_id is 1 for the BFS tree; if root is None or isolated, lists are minimal.
    """
    if root is None or root not in adj:
        # isolated or invalid root; return trivial
        return ([(root, None, 0.0, 0, 1)] if root else []), []

    # sort neighbors deterministically: by descending weight, then lex node id
    for u in list(adj.keys()):
        adj[u].sort(key=lambda t: (-(t[1]), (float('inf') if not math.isfinite(t[2]) else t[2]), t[0]))

    seen = {root}
    parent: Dict[str, Optional[str]] = {root: None}
    depth: Dict[str, int] = {root: 0}
    edge_w: Dict[Tuple[str, str], float] = {}
    edge_ev: Dict[Tuple[str, str], float] = {}
    edge_bs: Dict[Tuple[str, str], float] = {}

    q: deque[str] = deque([root])
    while q:
        u = q.popleft()
        for t in adj.get(u, []):
            v, w = t[0], t[1]
            ev = t[2] if len(t) > 2 else float('nan')
            bs = t[3] if len(t) > 3 else float('nan')
            if v in seen:
                continue
            seen.add(v)
            parent[v] = u
            depth[v] = depth[u] + 1
            edge_ev[(u, v)] = ev
            edge_bs[(u, v)] = bs
            edge_ev[(u, v)] = ev
            edge_bs[(u, v)] = bs
            edge_w[(u, v)] = w
            q.append(v)

    nodes_out: List[Tuple[str, Optional[str], float, int, int]] = []
    edges_out: List[Tuple[str, str, float, float, int]] = []

    comp_id = 1
    # root
    nodes_out.append((root, None, 0.0, 0, comp_id))

    # children
    for n, p in parent.items():
        if n == root:
            continue
        w = edge_w.get((p, n), edge_w.get((n, p), 0.0))
        nodes_out.append((n, p, w, depth[n], comp_id))
        edges_out.append((p, n, w, edge_ev.get((p, n), float('nan')), edge_bs.get((p, n), float('nan')), comp_id))

    # Sort outputs nicely
    nodes_out.sort(key=lambda x: (x[3], x[0]))  # by depth then node id
    edges_out.sort(key=lambda x: (depth.get(x[0], 0), -x[2], x[0], x[1]))  # parent depth, weight desc
    return nodes_out, edges_out



# === Added: BFS levels → Cytoscape Preset positions ===
from collections import defaultdict


# === Injected helpers: Cytoscape Preset positions + nodes TSV augmentation ===
from collections import defaultdict
import csv, gzip
from pathlib import Path as _Path_inj

def _compute_positions_from_nodes_out(nodes_out, x_gap=80.0, y_gap=120.0):
    """
    nodes_out: list of (node, parent, weight_to_parent, depth, component_id)
    Returns: dict node -> (x, y) layered by BFS depth.
    """
    levels = defaultdict(list)
    for node, parent, wtp, depth, comp_id in nodes_out:
        try:
            d = int(depth)
        except Exception:
            d = 0
        levels[d].append((node, parent))

    max_level = max(levels.keys()) if levels else 0

    # Deterministic order within each level: by parent then node id
    order_per_level = {}
    for d in range(max_level + 1):
        row = levels.get(d, [])
        row.sort(key=lambda t: ("" if t[1] is None else str(t[1]), str(t[0])))
        order_per_level[d] = [n for n, _ in row]

    pos = {}
    for d in range(max_level + 1):
        row = order_per_level.get(d, [])
        offset = -(len(row) - 1) / 2.0 if row else 0.0
        for i, n in enumerate(row):
            x = (offset + i) * x_gap
            y = -d * y_gap
            pos[n] = (x, y)
    return pos

def _augment_nodes_tsv_with_positions(nodes_path, positions):
    """
    In-place add 'X Location' and 'Y Location' columns to a nodes TSV/TSV.GZ.
    'positions' is a dict: node_id -> (x, y).
    Detects node id column among: 'node', 'name', 'uniprot_id'.
    """
    nodes_path = _Path_inj(nodes_path)
    is_gz = str(nodes_path).endswith(".gz")
    opener = gzip.open if is_gz else open

    with opener(nodes_path, "rt", encoding="utf-8", newline="") as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header in {nodes_path}")
        # Guess node id column
        if "node" in reader.fieldnames:
            id_col = "node"
        elif "name" in reader.fieldnames:
            id_col = "name"
        elif "uniprot_id" in reader.fieldnames:
            id_col = "uniprot_id"
        else:
            raise ValueError(f"Could not find a node id column in header: {reader.fieldnames}")

        new_fields = list(reader.fieldnames)
        if "X Location" not in new_fields:
            new_fields.append("X Location")
        if "Y Location" not in new_fields:
            new_fields.append("Y Location")

        rows = []
        for row in reader:
            nid = row.get(id_col, "")
            xy = positions.get(nid)
            if xy:
                row["X Location"] = xy[0]
                row["Y Location"] = xy[1]
            else:
                row["X Location"] = ""
                row["Y Location"] = ""
            rows.append(row)

    with opener(nodes_path, "wt", encoding="utf-8", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=new_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
# === End injected helpers ===


def _compute_positions_from_nodes_out(
    nodes_out: list[tuple[str, object, float, int, int]],
    x_gap: float = 80.0,
    y_gap: float = 120.0,
) -> dict[str, tuple[float, float]]:
    """
    Given nodes_out = [(node, parent, weight_to_parent, depth, component_id), ...],
    compute layered coordinates suitable for Cytoscape Preset layout.
    """
    # group by depth (level)
    levels = defaultdict(list)
    for node, parent, wtp, depth, comp_id in nodes_out:
        levels[int(depth)].append((node, parent))

    max_level = max(levels.keys()) if levels else 0

    # simple sibling ordering: sort by parent then node id for determinism
    order_per_level = {}
    for d in range(max_level + 1):
        row = levels.get(d, [])
        row.sort(key=lambda t: (str(t[1]) if t[1] is not None else "", str(t[0])))
        order_per_level[d] = [n for n, p in row]

    # assign coordinates centered per level
    pos = {}
    for d in range(max_level + 1):
        row = order_per_level.get(d, [])
        offset = -(len(row) - 1) / 2.0 if row else 0.0
        for i, n in enumerate(row):
            x = (offset + i) * x_gap
            y = -d * y_gap
            pos[n] = (x, y)
    return pos

def _write_cytoscape_positions_csv(
    pos: dict[str, tuple[float, float]],
    path: "Path",
    y_gap: float = 120.0,
) -> None:
    """
    Write CSV with columns Cytoscape Preset understands:
      name, X Location, Y Location, level
    """
    import csv
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["name", "X Location", "Y Location", "level"])
        for n, (x, y) in pos.items():
            level = int(round(-y / y_gap))
            w.writerow([n, x, y, level])

def main():
    nodes_path = NETWORK_DIR / "nodes.tsv.gz"
    edges_path = NETWORK_DIR / "edges_unique.tsv.gz"
    if not nodes_path.exists():
        raise SystemExit(f"Missing: {nodes_path}")
    if not edges_path.exists():
        raise SystemExit(f"Missing: {edges_path}")

    _ensure_outdir(OUT_DIR, overwrite=OVERWRITE_OUTPUT)

    # 1) Load node attributes and thresholded graph
    print("1) Load node attributes and thresholded graph")
    attrs = _read_nodes(nodes_path)
    adj   = _read_weighted_edges(edges_path, IDENTITY_THRESHOLD, EVALUE_MAX)

    # 2) Choose start node (prefer START_SUPERKINGDOM; else fallback to allowed SKs)
    print("2) Choose start node (prefer START_SUPERKINGDOM; else fallback to allowed SKs)")
    start = _choose_start_node(adj, attrs)
    if start is None:
        raise SystemExit(
            "No eligible start node found. Either no nodes match the chosen/allowed superkingdoms "
            f"(chosen={START_SUPERKINGDOM!r}, allowed={sorted(ALLOWED_SUPERKINGDOMS)}) "
            f"or no edges meet the thresholds (identity >= {IDENTITY_THRESHOLD}% and evalue <= {EVALUE_MAX})."
        )

    # 3) BFS go-forward traversal from 'start' (acyclic, no revisits)
    print("3) BFS go-forward traversal from 'start' (acyclic, no revisits)")
    nodes_out, edges_out = _bfs_go_forward(adj, start)

    # 4) Write outputs
    print("4) Write outputs")
    nodes_out_path = OUT_DIR / (f"go_forward_nodes_" + _format_run_suffix(IDENTITY_THRESHOLD, EVALUE_MAX, RUN_SUFFIX, FILTER_MODE, BITSCORE_MIN) + ".tsv.gz")
    edges_out_path = OUT_DIR / (f"go_forward_edges_" + _format_run_suffix(IDENTITY_THRESHOLD, EVALUE_MAX, RUN_SUFFIX, FILTER_MODE, BITSCORE_MIN) + ".tsv.gz")
    summary_path   = OUT_DIR / (f"go_forward_summary_" + _format_run_suffix(IDENTITY_THRESHOLD, EVALUE_MAX, RUN_SUFFIX, FILTER_MODE, BITSCORE_MIN) + ".txt")

    with gzip.open(nodes_out_path, "wt", encoding="utf-8") as fh:
        fh.write("node\tparent\tweight_to_parent\tdepth\tcomponent_id\tsuperkingdom\tgene\torganism\trank\tdegree_unique\n")
        for node, parent, wt, depth, comp in nodes_out:
            a = attrs.get(node)
            sk = a.superkingdom if a else ""
            gene = a.gene if a else ""
            org = a.organism if a else ""
            rank = a.rank if a else ""
            degu = a.degree_unique if a else 0
            wt_str = "" if parent is None else f"{wt:.3f}"
            fh.write(f"{node}\t{parent or ''}\t{wt_str}\t{depth}\t{comp}\t{sk}\t{gene}\t{org}\t{rank}\t{degu}\n")

    with gzip.open(edges_out_path, "wt", encoding="utf-8") as fh:
        fh.write("parent\tchild\tweight\tevalue\tbitscore\tcomponent_id\n")
        for p, c, w, ev, bs, comp in edges_out:
            fh.write(f"{p}\t{c}\t{w:.3f}\t{(ev if math.isfinite(ev) else float('inf'))}\t{(bs if math.isfinite(bs) else '')}\t{comp}\n")

    # 5) Summary
    start_attr = attrs.get(start)
    start_sk = start_attr.superkingdom if start_attr else ""
    with open(summary_path, "w", encoding="utf-8") as sf:
        sf.write("=== Go-forward network summary ===\n")
        sf.write(f"Edges file          : {edges_path}\n")
        sf.write(f"Nodes file          : {nodes_path}\n")
        sf.write(f"Identity threshold  : {IDENTITY_THRESHOLD:.3f}%\n")
        sf.write(f"Max e-value         : {EVALUE_MAX}\n")
        sf.write(f"Preferred start SK  : {START_SUPERKINGDOM}\n")
        sf.write(f"Allowed SK labels   : {', '.join(sorted(ALLOWED_SUPERKINGDOMS))}\n")
        sf.write("\n")
        sf.write(f"Start node          : {start}\n")
        sf.write(f"Start superkingdom  : {start_sk}\n")
        sf.write(f"Thresholded degree  : {len(adj.get(start, []))}\n")
        sf.write("\n")
        sf.write(f"Nodes in tree       : {len(nodes_out)}\n")
        sf.write(f"Edges in tree       : {len(edges_out)}\n")
        sf.write(f"Outputs:\n  {nodes_out_path}\n  {edges_out_path}\n")

    
    # 5) (Added) Add X/Y positions into the nodes TSV in place
    try:
        positions = _compute_positions_from_nodes_out(nodes_out, x_gap=80.0, y_gap=120.0)
        _augment_nodes_tsv_with_positions(nodes_out_path, positions)
        print(" [info] Added X/Y positions to nodes file:", nodes_out_path)
    except Exception as _pos_err:
        print("[warn] Failed to add positions to nodes file:", _pos_err)
    print("[done]")
    print(" Start node:", start, f"({_norm_sk(start_sk)})")
    print(" Wrote:")
    print("  ", nodes_out_path)
    print("  ", edges_out_path)
    print("  ", summary_path)


if __name__ == "__main__":
    main()
