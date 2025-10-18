# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3-listing_uniprot_merge_v5_gz_both_fast_imports.py
r16-2025-09-13

- Imports your parsers:
    from listing_parse import parse_ls_uniprot_scores
    from uniprot_parse import parse_uniprot_entry
- Streams BOTH Swiss-Prot and TrEMBL (.dat or .dat.gz) without full decompression
- Uses pigz for multi-core gzip streaming (if installed), else falls back to gzip.open
- Fast AC-line prefilter, then full parse via your uniprot_parse
- TSV header = preferred prefix + union of all other keys in first-seen order
- Ensures:
    * 'type' column comes from your listing (even when listing had a secondary accession)
    * 'protein' = recommended name only
    * 'gene'    = primary gene only
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple
import io
import gzip
import os
import re
import sys
import json
import shutil
import subprocess
from datetime import datetime

# External parsers (YOUR modules)
from listing_parse import parse_ls_uniprot_scores
from uniprot_parse import parse_uniprot_entry

try:
    from tqdm import tqdm
    _HAS_TQDM = True
except Exception:
    _HAS_TQDM = False


# --------------------
# Config / behavior toggles
# --------------------
# Dedup behavior across Swiss-Prot and TrEMBL:
#   None            -> keep both if the same accession appears in both
#   "prefer_swiss"  -> drop TrEMBL when accession already emitted from Swiss-Prot
#   "prefer_trembl" -> drop Swiss-Prot when accession already emitted from TrEMBL
DEDUP_MODE: Optional[str] = None

# TSV header always starts with this list (in this exact order).
PREFERRED_HEADER_PREFIX = [
    "rank", "tm", "coverage", "type",
    "protein",
    "gene",
    "protein_existence",
    "superkingdom",
    "organism",
    "lineage",
    "primary_accession",
    "interpro_ids",
    "source",
    "status",
    "id_line",
    "sequence_aa"
]


# --------------------
# Logging helpers
# --------------------
def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def info(msg: str) -> None:
    print(f"[{ts()}] [INFO] {msg}", flush=True)


# --------------------
# Path + opener helpers (pigz-enabled)
# --------------------
def resolve_dat_or_gz(pathlike: Path | str) -> Path:
    """
    If the .dat path doesn’t exist but a .dat.gz exists, return the .gz path; else return original.
    """
    p = Path(pathlike)
    if p.exists():
        return p
    s = str(p)
    if not s.endswith(".gz"):
        gz = Path(s + ".gz")
        if gz.exists():
            return gz
    return p

_PIGZ_DEBUG_SHOWN = False # print once if DEBUG_PIGZ=1

def open_text_auto(pathlike: Path) -> io.TextIOBase:
    """
    Return a text-mode file handle for a path that may be gzip-compressed (.gz) or plain.
    If the path ends with .gz and 'pigz' is available, use pigz -dc -p N for multi-core streaming.
    Env (optional):
      PIGZ_BIN      -> full path to pigz binary
      PIGZ_THREADS  -> integer threads; default ~half of CPU cores
      DEBUG_PIGZ    -> "1" to print which path is used (once)
    """
    global _PIGZ_DEBUG_SHOWN

    def _pigz_probe():
        import os, shutil
        print("\n[probe] DEBUG_PIGZ =", os.environ.get("DEBUG_PIGZ"))
        print("[probe] PIGZ_BIN    =", os.environ.get("PIGZ_BIN"))
        print("[probe] which pigz  =", shutil.which("pigz"))
    _pigz_probe()

    def _dbg(msg: str):
        global _PIGZ_DEBUG_SHOWN
        val = (os.environ.get("DEBUG_PIGZ") or "").lower()
        if val in {"1", "true", "yes", "all"}:
            # print once per process unless DEBUG_PIGZ=all
            if val == "all" or not _PIGZ_DEBUG_SHOWN:
                print(f"[pigz-debug pid={os.getpid()}] {msg}", file=sys.stderr, flush=True)
                _PIGZ_DEBUG_SHOWN = True

    p = Path(pathlike)
    s = str(p)
    if not s.endswith(".gz"):
        _dbg("using plain open() (not a .gz)")
        return p.open("r", encoding="utf-8", newline="")

    # Try pigz
    pigz_bin = os.environ.get("PIGZ_BIN") or shutil.which("pigz")
    if pigz_bin:
        n = os.environ.get("PIGZ_THREADS")
        if n is None:
            try:
                import multiprocessing as mp
                cores = mp.cpu_count() or 2
            except Exception:
                cores = 2
            n = str(max(2, cores // 2))  # default ~half cores to leave CPU for parsing
        try:
            proc = subprocess.Popen(
                [pigz_bin, "-dc", "-p", n, s],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
            )
            if proc.stdout is None:
                raise RuntimeError("pigz stdout is None")
            _dbg(f'using pigz at "{pigz_bin}" with {n} threads')
            return io.TextIOWrapper(proc.stdout, encoding="utf-8", newline="")
        except Exception as e:
            _dbg(f"pigz failed ({e!r}); falling back to gzip.open()")

    # Fallback: built-in gzip (single-core)
    _dbg("using builtin gzip.open()")
    return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", newline="")


# --------------------
# Entry streaming + fast prefilter
# --------------------
# Matches:
#  - legacy 6-char: [OPQ] + digit + 3 alnum + digit  OR  [A-NR-Z] + digit + 3 alnum + digit
#  - modern 10-char: A0 + 8 alnum  (covers A0A… family used by UniProtKB)
_UNIPROT_AC_REGEX = re.compile(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|A0[A-Z0-9]{8})\b")

def iter_uniprot_entries(path: Path | str) -> Iterator[List[str]]:
    """
    Yield lists of lines for each UniProtKB flat-file entry from `path`.
    Entry delimiter: a line that starts with '//'
    Works for .dat and .dat.gz (pigz-enabled via open_text_auto).
    """
    with open_text_auto(resolve_dat_or_gz(path)) as fh:
        buf: List[str] = []
        for line in fh:
            if line.startswith("//"):
                if buf:
                    yield buf
                    buf = []
            else:
                buf.append(line)
        if buf:
            yield buf

def find_first_matched_accession(entry_lines: List[str], wanted_set: Set[str]) -> Optional[str]:
    """
    Return the FIRST accession in AC lines that is present in wanted_set, else None.
    This lets us merge listing metadata by the *matched* accession (which could be secondary),
    ensuring fields like 'type' are not lost.
    """
    for ln in entry_lines:
        if ln.startswith("AC "):
            for m in _UNIPROT_AC_REGEX.finditer(ln):
                ac = m.group(0)
                if ac in wanted_set:
                    return ac
    return None



# --------------------
# TSV writer with dynamic header + SORT (rank asc, tm desc, coverage desc)
# --------------------
import tempfile
import heapq
import math
import os

def _to_int_safe(v, default=10**9):
    try:
        if v is None or v == "":
            return default
        return int(str(v).strip())
    except Exception:
        return default

def _to_float_safe(v, default=-math.inf):
    try:
        if v is None or v == "":
            return default
        return float(str(v).strip())
    except Exception:
        return default

def _sort_key(row: dict):
    # rank asc, tm desc, coverage desc
    r = _to_int_safe(row.get("rank"), default=10**9)
    tm = _to_float_safe(row.get("tm"), default=-float("inf"))
    cov = _to_float_safe(row.get("coverage"), default=-float("inf"))
    return (r, -tm, -cov)

def write_dynamic_tsv_sorted(out_path: Path, rows: Iterable[Dict[str, str]], preferred_prefix: List[str]) -> None:
    """
    Streams rows to JSONL to learn keys, then sorts rows externally.
    Header = preferred_prefix + union(other keys) in first-seen order.
    Prefix order always overrides (no dupes).
    """
    tmp_dir = out_path.parent
    stage_jsonl = out_path.with_suffix(out_path.suffix + ".stage.jsonl")
    chunk_files = []

    # Learn keys (skipping any in prefix)
    prefix_set = set(preferred_prefix)
    seen_nonpref: set[str] = set()
    nonpref_order: list[str] = []

    with stage_jsonl.open("w", encoding="utf-8", newline="") as op:
        for row in rows:
            # ensure preferred columns exist
            for k in preferred_prefix:
                row.setdefault(k, "")
            # learn non-preferred keys only
            for k in row.keys():
                if k in prefix_set:
                    continue
                if k not in seen_nonpref:
                    seen_nonpref.add(k)
                    nonpref_order.append(k)
            op.write(json.dumps(row, ensure_ascii=False) + "\n")

    # Final header with strict prefix-first rule
    seen = set()
    header = []
    for k in preferred_prefix:
        if k not in seen:
            header.append(k); seen.add(k)
    for k in nonpref_order:
        if k not in prefix_set and k not in seen:
            header.append(k); seen.add(k)

    # --- existing chunk sort & merge code below this line stays the same ---
    CHUNK_ROWS = int(os.environ.get("SORT_CHUNK_ROWS", "250000"))
    buf = []

    def _sort_key(row: dict):
        # rank asc, tm desc, coverage desc
        r = _to_int_safe(row.get("rank"), default=10**9)
        tm = _to_float_safe(row.get("tm"), default=-float("inf"))
        cov = _to_float_safe(row.get("coverage"), default=-float("inf"))
        return (r, -tm, -cov)


    def flush_chunk(buf_rows, idx):
        buf_rows.sort(key=_sort_key)
        cf = Path(tempfile.mkstemp(prefix=f"uniprot_chunk_{idx:04d}_", suffix=".jsonl", dir=tmp_dir)[1])
        with cf.open("w", encoding="utf-8", newline="") as cfp:
            for r in buf_rows:
                cfp.write(json.dumps(r, ensure_ascii=False) + "\n")
        return cf

    with stage_jsonl.open("r", encoding="utf-8") as ip:
        for line in ip:
            buf.append(json.loads(line))
            if len(buf) >= CHUNK_ROWS:
                chunk_files.append(flush_chunk(buf, len(chunk_files)))
                buf = []
        if buf:
            chunk_files.append(flush_chunk(buf, len(chunk_files)))
            buf = []

    def iter_chunk(cf_path: Path):
        with cf_path.open("r", encoding="utf-8") as fh:
            for ln in fh:
                yield json.loads(ln)

    iters = [iter_chunk(cf) for cf in chunk_files]
    heap = []
    for i, it in enumerate(iters):
        try:
            r = next(it)
            heap.append((_sort_key(r), i, r))
        except StopIteration:
            pass
    heapq.heapify(heap)

    with out_path.open("w", encoding="utf-8", newline="") as tsv:
        tsv.write("\t".join(header) + "\n")
        while heap:
            _, src_i, row = heapq.heappop(heap)
            tsv.write("\t".join("" if row.get(k) is None else str(row.get(k)) for k in header) + "\n")
            try:
                nxt = next(iters[src_i])
                heapq.heappush(heap, (_sort_key(nxt), src_i, nxt))
            except StopIteration:
                pass

    try:
        stage_jsonl.unlink()
    except Exception:
        pass
    for cf in chunk_files:
        try:
            cf.unlink()
        except Exception:
            pass

# --------------------
# Helpers to enforce protein/gene rules on the parsed dict
# --------------------
def _extract_recommended_protein(meta: Dict[str, object]) -> str:
    # Common key patterns that might hold the recommended name
    for k in ("recommended_name", "protein_recommended", "protein_recommended_name", "protein_name_recommended"):
        v = meta.get(k)
        if v:
            return str(v)
    # Fallbacks
    v = meta.get("protein")
    if v:
        # If parser put multiple names in one string, pick the first token reasonably
        s = str(v)
        for sep in [" | ", ";", " (", " [", " {", "\t"]:
            if sep in s:
                return s.split(sep)[0].strip()
        return s.strip()
    return ""

def _extract_primary_gene(meta: Dict[str, object]) -> str:
    for k in ("primary_gene", "gene_primary", "gene_name_primary"):
        v = meta.get(k)
        if v:
            return str(v)
    # Fallbacks
    v = meta.get("gene")
    if v:
        s = str(v)
        # Split on common separators; keep first token
        for sep in [",", ";", " "]:
            if sep in s:
                return s.split(sep)[0].strip()
        return s.strip()
    return ""


# --------------------
# Single-file scanner (uses your parsers)
# --------------------
def scan_uniprot_file(
    uni_path: Path,
    wanted: Set[str],
    listing_meta: Dict[str, Dict[str, str]],
    source_label: str,
    columns: List[str],  # kept for back-compat; dynamic writer ignores this
    dedup_seen: Optional[Set[str]] = None,
    progress_desc: str = "Entries scanned"
) -> Iterator[Dict[str, str]]:
    """
    Stream a UniProtKB .dat(.gz) file and yield *dict* rows for matches.
    Row keys are fields from parse_uniprot_entry() plus:
      - any listing fields for the *matched accession* (may be secondary!)
      - 'source' -> 'Swiss-Prot' or 'TrEMBL'
    If dedup_seen is provided, drop entries whose primary_accession was already emitted.
    """
    pbar = tqdm(disable=not _HAS_TQDM, desc=progress_desc, unit="entry")
    try:
        for entry_lines in iter_uniprot_entries(uni_path):
            if _HAS_TQDM:
                pbar.update(1)

            matched_id = find_first_matched_accession(entry_lines, wanted)
            if not matched_id:
                continue

            # Full parse via your external parser
            entry_text = ''.join(entry_lines)
            meta_map = parse_uniprot_entry(entry_text)
            if not meta_map:
                continue

            # Expect a single primary accession key
            primary_ac, meta = next(iter(meta_map.items()))
            if not primary_ac:
                continue

            # Dedup if requested
            if dedup_seen is not None:
                if primary_ac in dedup_seen:
                    continue
                dedup_seen.add(primary_ac)

            # Start row with listing fields for the *matched accession* (so 'type' is kept)
            row: Dict[str, str] = {}
            lm = listing_meta.get(matched_id, {})
            if lm:
                row.update({k: ("" if v is None else str(v)) for k, v in lm.items()})
            # Make sure 'type' is present even if lm was empty
            row.setdefault("type", lm.get("type", ""))

            # Identify data source
            row["source"] = source_label

            # Add parsed fields (don’t clobber same-named listing fields)
            for k, v in meta.items():
                row.setdefault(k, "" if v is None else str(v))

            # Ensure primary_accession present
            row.setdefault("primary_accession", primary_ac)

            # Enforce: protein → recommended only; gene → primary only
            prot = _extract_recommended_protein(meta)
            if prot:
                row["protein"] = prot
            gene = _extract_primary_gene(meta)
            if gene:
                row["gene"] = gene

            yield row
    finally:
        if _HAS_TQDM:
            pbar.close()


# --------------------
# Pipeline (no CLI; edit paths in __main__)
# --------------------
def run_uniprot_merge_pipeline(
    listing_path: Path,
    sprot_path: Path,
    trembl_path: Path,
    out_tsv_path: Path,
    columns: List[str] = PREFERRED_HEADER_PREFIX,  # ignored by dynamic writer
    parse_listing_mode: str = "filename_tokens",   # ignored here; your parser decides
    limit: int = 0,  # 0 = no cap on IDs; else cap wanted IDs (debug)
) -> None:
    # Step 1: build the wanted set (and listing metadata) using your module
    info("Step 1/4: Parsing listing with external parser ...")
    listing_meta = parse_ls_uniprot_scores(listing_path)
    wanted = set(listing_meta.keys())
    if limit and limit > 0:
        wanted = set(list(sorted(wanted))[:limit])  # deterministic cap
    info(f"Parsed {len(wanted):,} UniProt IDs from listing.")

    # Optional dedup set depending on preference
    dedup_seen: Optional[Set[str]] = set() if DEDUP_MODE in ("prefer_swiss", "prefer_trembl") else None

    # Step 2: scan Swiss-Prot
    sprot_path = resolve_dat_or_gz(sprot_path)
    info(f"Step 2/4: Scanning Swiss-Prot file: {sprot_path}")
    swiss_iter = scan_uniprot_file(
        sprot_path,
        wanted,
        listing_meta,
        source_label="Swiss-Prot",
        columns=columns,
        dedup_seen=(dedup_seen if DEDUP_MODE == "prefer_swiss" else None),
        progress_desc="Entries scanned (Swiss-Prot)"
    )

    # Step 2b: scan TrEMBL
    trembl_path = resolve_dat_or_gz(trembl_path)
    info(f"Step 2b/4: Scanning TrEMBL file: {trembl_path}")
    trembl_iter = scan_uniprot_file(
        trembl_path,
        wanted,
        listing_meta,
        source_label="TrEMBL",
        columns=columns,
        dedup_seen=(dedup_seen if DEDUP_MODE == "prefer_trembl" else None),
        progress_desc="Entries scanned (TrEMBL)"
    )

    # Step 3: (optional) ranking — streaming keeps order; add sorting if you need global order.
    info("Step 3/4: Ranking rows (if applicable) ...")

    # Step 4: write TSV with dynamic header (preferred prefix + union of others)
    info(f"Step 4/4: Writing TSV -> {out_tsv_path}")

    def dict_row_gen():
        if DEDUP_MODE == "prefer_trembl":
            for r in trembl_iter:
                yield r
            for r in swiss_iter:
                yield r
        else:
            for r in swiss_iter:
                yield r
            for r in trembl_iter:
                yield r

    write_dynamic_tsv_sorted(out_tsv_path, dict_row_gen(), PREFERRED_HEADER_PREFIX)
    info("Done.")


# --------------------
# Main (edit paths here; no CLI)
# --------------------
if __name__ == "__main__":
    # EDIT THESE PATHS for your environment:
    LISTING_PATH = Path("/Users/danisom/Desktop/local_informatics/2023.11.12.7tmps_results-filtered_hits_listing.txt")

    SPROT_PATH  = Path("/Volumes/isomlab_1/bulk_uniprot/uniprot_sprot.dat.gz")    # or .dat
    TREMBL_PATH = Path("/Volumes/isomlab_1/bulk_uniprot/uniprot_trembl.dat.gz")   # or .dat

    OUTPUT_TSV  = Path("/Volumes/isomlab_1/bulk_uniprot/filtered_merged.tsv")

    # Optional: help pigz pick sensible threads (or set in shell: export PIGZ_THREADS=12)
    os.environ.setdefault("PIGZ_BIN", "/opt/anaconda3/envs/tqdm/bin/pigz")
    os.environ.setdefault("PIGZ_THREADS", "12")  # e.g., good start for 14-core M4 Max

    # Optional: pigz debug once at startup (prints which path is used)
    os.environ.setdefault("DEBUG_PIGZ", "1")

    try:
        run_uniprot_merge_pipeline(
            listing_path=LISTING_PATH,
            sprot_path=SPROT_PATH,
            trembl_path=TREMBL_PATH,
            out_tsv_path=OUTPUT_TSV,
            columns=PREFERRED_HEADER_PREFIX,  # kept for signature; dynamic writer ignores beyond prefix
            parse_listing_mode="filename_tokens",
            limit=0,
        )
    except FileNotFoundError as e:
        sys.stderr.write(f"[ERROR] {e}\n")
        sys.exit(1)
    except KeyboardInterrupt:
        sys.stderr.write("[ERROR] Interrupted by user.\n")
        sys.exit(2)
    except Exception as e:
        sys.stderr.write(f"[ERROR] {type(e).__name__}: {e}\n")
        raise
