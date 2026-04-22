# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reduce_for_cytoscape_strict.py
------------------------------
Edit the CONFIG block and run:  python reduce_for_cytoscape_strict.py
(No command-line flags.)

Purpose
- Reduce a very large network to a Cytoscape-friendly subset while keeping the most
  informative structure for figures: connectivity + strongest relationships.
- Works on Windows/macOS/Linux with stock Python 3 (no extra packages).

Key ideas
- Keep only the largest components until a global node budget.
- For each kept component:
  * Build a Maximum Spanning Forest (MSF) using edge weights (strongest backbone).
  * If the component is still too big, prune MSF leaves by node "importance" score
    until a per-component cap is met. Score = degree_unique (if present) + sum of
    top-3 incident edge weights (local density proxy).
  * Optionally add a few (K) strongest non-backbone edges per node to restore local density.
- Optional: drop weak edges early (MIN_WEIGHT) and low-degree nodes at the end (MIN_DEGREE).
- Writes two TSVs for Cytoscape; if nodes have 'depth', also writes x,y for Preset layout.

Inputs (auto-detected)
- Nodes TSV (.tsv or .tsv.gz): must contain 'uid' OR 'node', and 'component_id'.
- Edges TSV (.tsv or .tsv.gz): either (src_uid, dst_uid, pident_max) or (parent, child, weight).

Cytoscape import
1) File → Import → Network from File… → <OUT_PREFIX>_edges.tsv (Source=source, Target=target, set Directed if you want arrows)
2) File → Import → Table from File… → <OUT_PREFIX>_nodes.tsv (To: Node Table; map 'uid' or 'node' → 'name')
3) Layout → Preset (X=x, Y=y) if present → Fit Content; then style and export SVG/PDF.
"""

import csv, gzip, heapq
from pathlib import Path
from collections import defaultdict, Counter, namedtuple

# ------------------------ CONFIG (EDIT THESE) ------------------------
NODES_PATH         = "/Users/danisom/go_forward_out_id30_b30_abe_cu/go_forward_nodes_id30_b30_abe_cu.tsv"      # .tsv or .tsv.gz; must have 'uid' or 'node' and 'component_id'
OUT_DIR            = Path(NODES_PATH).parent  # write outputs next to nodes file
EDGES_PATH         = "/Users/danisom/go_forward_out_id30_b30_abe_cu/go_forward_edges_id30_b30_abe_cu.tsv"      # .tsv or .tsv.gz; (src_uid,dst_uid,pident_max) or (parent,child,weight)
OUT_PREFIX         = "cyto_reduced_strict"


TARGET_MAX_NODES   = 12000   # total node budget across largest components (12000 is default)
MAX_PER_COMPONENT  = 2000    # hard cap per component (enforced by pruning) (8000 is default)
MIN_WEIGHT         = 0.0     # drop edges below this before reduction (e.g., 30.0)
PER_NODE_EXTRA     = 1       # add up to K strongest extra edges per node after MSF (0..2 typical)
MIN_DEGREE         = 0       # drop nodes with final degree < this (e.g., 2 for cleaner figs)

# Auto suffix including key config parameters for traceable outputs
OUT_SUFFIX = f"_tm{int(TARGET_MAX_NODES)}_mc{int(MAX_PER_COMPONENT)}_mw{MIN_WEIGHT:g}_pe{int(PER_NODE_EXTRA)}_md{int(MIN_DEGREE)}"
# --------------------------------------------------------------------

Edge = namedtuple("Edge", "u v w ev raw_dir")  # raw_dir: 1 if original row was directed u->v, 0 if undirected source

def gz_or_tsv_reader(path: Path):
    # Open and build a DictReader, but normalize headers and row keys
    if str(path).lower().endswith(".gz"):
        fh = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        fh = open(path, "r", encoding="utf-8", newline="")
    rdr = csv.DictReader(fh, delimiter="\t")
    # Normalize header names: strip spaces and BOM
    if rdr.fieldnames:
        normalized = [ (h or "").lstrip("\ufeff").strip() for h in rdr.fieldnames ]
        rdr.fieldnames = normalized
    # Wrap the iterator to normalize each row's keys as well
    class _RowNormalizer:
        def __init__(self, reader):
            self._reader = reader
            self.fieldnames = reader.fieldnames or []
        def __iter__(self):
            for row in self._reader:
                yield { (k or "").lstrip("\ufeff").strip(): v for k, v in row.items() }
    return fh, _RowNormalizer(rdr)

def detect_node_id_col(header):
    # normalize once
    norm = { (h or "").lstrip("\ufeff").strip().lower(): h for h in header }
    # prefer canonical names
    for cand in ("uid","node","id","node_id","uid_id"):
        if cand in norm:
            return norm[cand]
    raise SystemExit("Nodes file must contain a node identifier column (uid/node/id).")

def detect_edge_cols(header):
    if all(k in header for k in ("src_uid","dst_uid","pident_max")):
        return ("src_uid","dst_uid","pident_max", False)
    if all(k in header for k in ("parent","child","weight")):
        return ("parent","child","weight", True)
    if all(k in header for k in ("src_uid","dst_uid")):
        return ("src_uid","dst_uid", None, False)
    if all(k in header for k in ("parent","child")):
        return ("parent","child", None, True)
    raise SystemExit("Edges file must contain (src_uid,dst_uid,pident_max) or (parent,child,weight).")

class UnionFind:
    def __init__(self, nodes=None):
        self.parent = {}
        self.rank = {}
        if nodes:
            for x in nodes:
                self.parent[x] = x
                self.rank[x] = 0
    def find(self, x):
        p = self.parent.get(x, x)
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent.get(x, x)
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb: return False
        if self.rank.get(ra,0) < self.rank.get(rb,0):
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank.get(ra,0) == self.rank.get(rb,0):
            self.rank[ra] = self.rank.get(ra,0) + 1
        return True

def pick_components(nodes_path: Path, target_max_nodes: int):
    fh, rdr = gz_or_tsv_reader(nodes_path)
    with fh:
        header = rdr.fieldnames or []
        nid_col = detect_node_id_col(header)
        if "component_id" not in header:
            raise SystemExit("Nodes file must contain 'component_id' column.")
        comp_by_node = {}
        nodes_by_comp = defaultdict(list)
        for row in rdr:
            nid = row[nid_col]
            comp = row["component_id"]
            comp_by_node[nid] = comp
            nodes_by_comp[comp].append(row)

    comps_sorted = sorted(nodes_by_comp.items(), key=lambda kv: len(kv[1]), reverse=True)
    kept, total = [], 0
    for comp, rows in comps_sorted:
        if kept and total >= target_max_nodes:
            break
        kept.append(comp)
        total += len(rows)
    print(f"[pick] Keeping {len(kept)} components (≈{total} nodes) within target={target_max_nodes}")
    return nid_col, comp_by_node, nodes_by_comp, kept  # keep order

from collections import namedtuple
Edge = namedtuple('Edge', ['u','v','w','ev','bs','raw_dir'])

def read_edges_grouped(edges_path: Path, comp_of, kept_components, min_weight: float):
    fh, rdr = gz_or_tsv_reader(edges_path)
    with fh:
        header = rdr.fieldnames or []
        s_col, t_col, w_col, is_directed = detect_edge_cols(header)
        ev_col = 'evalue' if 'evalue' in header else None
        bs_col = 'bitscore' if 'bitscore' in header else None
        ev_col = 'evalue' if 'evalue' in header else None
        bs_col = 'bitscore' if 'bitscore' in header else None
        edges_by_comp = defaultdict(list)
        kept = 0
        for row in rdr:
            u = row[s_col]; v = row[t_col]
            cu = comp_of.get(u); cv = comp_of.get(v)
            if not cu or not cv or cu != cv or cu not in kept_components:
                continue
            w = 1.0
            if w_col:
                try:
                    w = float(row[w_col])
                except:
                    continue
                if w < min_weight:
                    continue
            raw_dir = 1 if is_directed else 0
            edges_by_comp[cu].append(Edge(
                u=u, v=v, w=w,
                ev=(
                    float('inf', bs=float('nan')) if (ev_col and ((row.get(ev_col) or '').strip().lower() in ('inf','infinity')))
                    else (float(row[ev_col]) if (ev_col and (row.get(ev_col) or '').strip()!='') else float('nan'))
                ),
                bs=(float(row[bs_col]) if (bs_col and (row.get(bs_col) or '').strip()!='') else float('nan')),
                raw_dir=raw_dir
            ))
            kept += 1
    print(f"[edges] Kept {kept:,} edges after component/weight filtering")
    return edges_by_comp

def maximum_spanning_forest(nodes_in_comp, edges):
    edges_sorted = sorted(edges, key=lambda e: e.w, reverse=True)
    uf = UnionFind(nodes_in_comp)
    adj = defaultdict(list)  # MSF adjacency
    picked = []
    for e in edges_sorted:
        a, b = e.u, e.v
        if uf.union(a, b):
            picked.append(e)
            adj[a].append((b, e.w))
            adj[b].append((a, e.w))
    for n in nodes_in_comp:
        adj[n]  # ensure key
    return adj, picked

def compute_node_scores(msf_adj, all_edges_in_comp, nodes_meta):
    # Sum top-3 incident weights (local density proxy)
    inc = defaultdict(list)
    for e in all_edges_in_comp:
        inc[e.u].append(e.w)
        inc[e.v].append(e.w)
    weighted = {n: sum(sorted(wts, reverse=True)[:3]) for n, wts in inc.items()}
    scores = {}
    for n in msf_adj.keys():
        du = 0.0
        meta = nodes_meta.get(n)
        if meta is not None:
            val = meta.get("degree_unique") if isinstance(meta, dict) else None
            if val is not None:
                try:
                    du = float(val)
                except:
                    du = 0.0
        scores[n] = du + weighted.get(n, 0.0)
    return scores

def prune_msf_to_cap(msf_adj, cap, scores):
    if cap <= 0:
        return set(), []
    # current degrees
    adj = {n: set(m for m,_ in nbrs) for n, nbrs in msf_adj.items()}
    deg = {n: len(nei) for n, nei in adj.items()}
    # min-heap of (score, node) for leaves
    heap = [(scores.get(n, 0.0), n) for n,d in deg.items() if d <= 1]
    heapq.heapify(heap)
    kept = set(adj.keys())

    def remove_leaf(x):
        kept.discard(x)
        for nb in list(adj[x]):
            adj[nb].discard(x)
            deg[nb] -= 1
            if deg[nb] == 1 and nb in kept:
                heapq.heappush(heap, (scores.get(nb, 0.0), nb))
        adj[x].clear()
        deg[x] = 0

    while len(kept) > cap and heap:
        _, node = heapq.heappop(heap)
        if node not in kept:
            continue
        remove_leaf(node)

    kept_edges = []
    for u in kept:
        for v in adj[u]:
            if u < v and v in kept:
                kept_edges.append((u, v))
    return kept, kept_edges

def add_top_extra_edges_per_node(all_edges_in_comp, kept_nodes, existing_pairs, per_node_extra):
    if per_node_extra <= 0:
        return []
    cand_by_node = defaultdict(list)
    for e in all_edges_in_comp:
        if e.u in kept_nodes and e.v in kept_nodes:
            a, b = (e.u, e.v) if e.u <= e.v else (e.v, e.u)
            cand_by_node[e.u].append((a,b,e))
            cand_by_node[e.v].append((a,b,e))
    added_pairs = set()
    added = []
    quota = Counter()
    for node, lst in cand_by_node.items():
        lst.sort(key=lambda t: t[2].w, reverse=True)
        taken = 0
        for a,b,e in lst:
            if quota[node] >= per_node_extra:
                break
            if (a,b) in existing_pairs or (a,b) in added_pairs:
                continue
            added_pairs.add((a,b))
            added.append(e)
            quota[node] += 1
            taken += 1
            if taken >= per_node_extra:
                break
    return added

def safe_int(x, default=0):
    try:
        return int(float(x))
    except:
        return default

def maybe_make_xy(nodes_rows):
    if not nodes_rows or "depth" not in nodes_rows[0].keys():
        return None
    by_depth = defaultdict(list)
    for r in nodes_rows:
        try:
            d = int(float(r.get("depth", "0") or 0))
        except:
            d = 0
        by_depth[d].append(r)
    X_SPACING = 120.0
    Y_SPACING = 18.0
    positions = {}
    for depth, rows in by_depth.items():
        rows_sorted = sorted(rows, key=lambda r: (-safe_int(r.get("degree_unique","0")), r.get("uid") or r.get("node")))
        h = len(rows_sorted)
        for idx, r in enumerate(rows_sorted):
            nid = r.get("uid") or r.get("node")
            y_index = idx - (h - 1) / 2.0
            x = depth * X_SPACING
            y = y_index * Y_SPACING
            positions[nid] = (x, y)
    return positions

def main():
    nodes_path = Path(NODES_PATH)
    edges_path = Path(EDGES_PATH)
    out_nodes = OUT_DIR / f"{OUT_PREFIX}{OUT_SUFFIX}_nodes.tsv"
    out_edges = OUT_DIR / f"{OUT_PREFIX}{OUT_SUFFIX}_edges.tsv"

    # [1] Nodes & components
    nid_col, comp_of, nodes_by_comp, kept_components_ordered = pick_components(nodes_path, TARGET_MAX_NODES)
    kept_components = set(kept_components_ordered)

    # build per-comp meta
    nodes_meta_by_comp = {}
    for comp, rows in nodes_by_comp.items():
        d = {}
        for r in rows:
            nid = r.get("uid") or r.get("node")
            d[nid] = r
        nodes_meta_by_comp[comp] = d

    # [2] Edges grouped
    edges_by_comp = read_edges_grouped(edges_path, comp_of, kept_components, MIN_WEIGHT)

    # [3] Per component: MSF → prune to cap → add extras
    final_edges = []
    kept_nodes_global = set()

    comp_sizes = {c: len(nodes_by_comp[c]) for c in kept_components}
    total_kept_nodes = sum(comp_sizes[c] for c in kept_components)

    for comp in kept_components_ordered:
        comp_edges = edges_by_comp.get(comp, [])
        if not comp_edges:
            continue
        nodes_in_comp = set()
        for e in comp_edges:
            nodes_in_comp.add(e.u); nodes_in_comp.add(e.v)

        msf_adj, msf_edges = maximum_spanning_forest(nodes_in_comp, comp_edges)

        proportional_cap = max(1, int(TARGET_MAX_NODES * (comp_sizes[comp] / max(1, total_kept_nodes))))
        cap = min(MAX_PER_COMPONENT, proportional_cap)

        scores = compute_node_scores(msf_adj, comp_edges, nodes_meta_by_comp.get(comp, {}))
        kept_nodes, kept_msf_pairs = prune_msf_to_cap(msf_adj, cap, scores)

        # record MSF edges
        for (a,b) in kept_msf_pairs:
            w_ab = max((e.w for e in comp_edges if (e.u==a and e.v==b) or (e.u==b and e.v==a)), default=1.0)
            ev_ab = min((e.ev for e in comp_edges if ((e.u==a and e.v==b) or (e.u==b and e.v==a)) and (e.ev==e.ev)), default=float('inf'))
            bs_ab = max((getattr(e, 'bs', float('nan')) for e in comp_edges if (e.u==a and e.v==b) or (e.u==b and e.v==a)), default=float('nan'))
            final_edges.append({
                "source": a, "target": b,
                "weight": w_ab, "evalue": ev_ab, "bitscore": bs_ab,
                "component_id": comp, "is_backbone": 1
            })

        # add extras among kept nodes
        existing_pairs = set((min(r["source"], r["target"]), max(r["source"], r["target"])) for r in final_edges if r["component_id"]==comp)
        added = add_top_extra_edges_per_node(comp_edges, kept_nodes, existing_pairs, PER_NODE_EXTRA)
        for e in added:
            a,b = (e.u, e.v) if e.u <= e.v else (e.v, e.u)
            if (a,b) not in existing_pairs:
                final_edges.append({
                    "source": e.u, "target": e.v, "weight": e.w,
                    "evalue": (getattr(e, 'ev', float('inf'))),
                    "bitscore": (getattr(e, 'bs', float('nan'))),
                    "component_id": comp, "is_backbone": 0
                })
                existing_pairs.add((a,b))
    kept_nodes_global |= kept_nodes

    if not final_edges:
        print("No edges survived filtering. Consider lowering MIN_WEIGHT or raising TARGET_MAX_NODES/MAX_PER_COMPONENT.")
        with out_nodes.open("w", encoding="utf-8", newline="") as fn, out_edges.open("w", encoding="utf-8", newline="") as fe:
            fn.write("uid\tcomponent_id\n")
            fe.write("source\ttarget\tweight\tcomponent_id\tis_backbone\n")
        return

    # [4] Optional final degree prune
    if MIN_DEGREE > 0:
        deg = Counter()
        for r in final_edges:
            deg[r["source"]] += 1; deg[r["target"]] += 1
        final_edges = [r for r in final_edges if deg[r["source"]] >= MIN_DEGREE and deg[r["target"]] >= MIN_DEGREE]
        kept_nodes_global = set()
        for r in final_edges:
            kept_nodes_global.add(r["source"]); kept_nodes_global.add(r["target"])

    # [5] Write edges TSV
    with out_edges.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["source","target","weight","evalue","bitscore","component_id","is_backbone"])
        for r in final_edges:
            w.writerow([r["source"], r["target"], r.get("weight",""), r.get("evalue",""), r.get("bitscore",""), r["component_id"], r["is_backbone"]])

    # [6] Write nodes TSV (preserve columns; add x,y if depth exists)
    kept_components_final = set(r["component_id"] for r in final_edges)

    # find a header
    hdr_nodes = None
    for comp in kept_components_ordered:
        if comp in kept_components_final and nodes_by_comp[comp]:
            hdr_nodes = list(nodes_by_comp[comp][0].keys())
            break
    if hdr_nodes is None:
        hdr_nodes = ["uid","component_id"]

    id_col = "uid" if "uid" in hdr_nodes else ("node" if "node" in hdr_nodes else detect_node_id_col(hdr_nodes))

    kept_rows = []
    for comp in kept_components_ordered:
        if comp not in kept_components_final:
            continue
        for r in nodes_by_comp[comp]:
            nid = r.get(id_col) or r.get("uid") or r.get("node")
            if nid in kept_nodes_global:
                kept_rows.append(r)

    positions = None
    if kept_rows and "depth" in kept_rows[0].keys():
        positions = maybe_make_xy(kept_rows)

    final_hdr = list(hdr_nodes)
    if "x" not in final_hdr: final_hdr.append("x")
    if "y" not in final_hdr: final_hdr.append("y")

    with out_nodes.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(final_hdr)
        for r in kept_rows:
            row_out = [r.get(col, "") for col in hdr_nodes]
            nid = r.get(id_col) or r.get("uid") or r.get("node")
            if positions and nid in positions:
                x, y = positions[nid]
            else:
                x, y = "", ""
            row_out += [x, y]
            w.writerow(row_out)

    print(f"Done. Wrote: {out_nodes} and {out_edges}")
    print("Cytoscape import:\n  Network: edges (source/target)\n  Table: nodes (map id→name)\n  Layout: Preset (x,y) if present; Fit Content.\n  Style: Color=superkingdom, Size=degree_unique, Edge width=weight.")

if __name__ == "__main__":
    main()
