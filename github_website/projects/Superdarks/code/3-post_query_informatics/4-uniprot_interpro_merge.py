# uniprot_interpro_merge_v11.py
# v11: Force UPPERCASE for 'protein', 'gene', and InterPro descriptions in output.
#
from pathlib import Path
from typing import Dict, List, Set, Tuple
import csv
import gzip
import io
import os
import re
import shutil
import subprocess
import sys

# ==========================
#  CONFIGURE YOUR FILE PATHS
# ==========================
INPUT_TSV = "/Volumes/isomlab_1/bulk_uniprot/filtered_merged.tsv"           # your source .tsv (like the excerpt)
PROTEIN2IPR_GZ = "/Volumes/isomlab_1/bulk_interpro/protein2ipr.dat.gz"      # InterPro protein2ipr flat file (gzipped)
OUTPUT_TSV = "/Volumes/isomlab_1/bulk_interpro/output_with_interpro.tsv"    # new file to write

# Optional caps for testing
TEST_LIMIT_GZ_LINES = 0               # 0 = disabled; stop after N gz lines (fast dry run)

# ==============
#  tqdm fallback
# ==============
try:
    from tqdm import tqdm
    def _tqdm(iterable=None, **kwargs):
        return tqdm(iterable=iterable, ncols=80, leave=False, **kwargs)
    TQDM = True
except Exception:
    TQDM = False
    def _tqdm(iterable=None, **kwargs):
        return iterable

def _info(msg: str):
    print(msg, file=sys.stderr, flush=True)

# ---------------------------
# pigz-backed streaming open
# ---------------------------
def open_text_auto(pathlike: Path, pigz_threads: int = 12) -> io.TextIOBase:
    p = Path(pathlike)
    s = str(p)
    if not s.endswith(".gz"):
        return p.open("r", encoding="utf-8", newline="")

    pigz_bin = os.environ.get("PIGZ_BIN") or shutil.which("pigz")
    if pigz_bin:
        threads = str(int(os.environ.get("PIGZ_THREADS", pigz_threads)))
        proc = subprocess.Popen(
            [pigz_bin, "-dc", "-p", threads, s],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            bufsize=1,
        )
        if proc.stdout is None:
            raise RuntimeError("pigz stdout is None")
        _info(f"[3/4] Using pigz at '{pigz_bin}' with {threads} threads")
        return proc.stdout
    else:
        _info("[3/4] pigz not found on PATH; falling back to gzip.open (single-core)")
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", newline="")

def count_data_rows(tsv_path: Path) -> int:
    total = 0
    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        _ = f.readline()  # header
        for _ in f:
            total += 1
    return total

def parse_tsv_collect(tsv_path: Path) -> Tuple[List[Dict[str,str]], Set[str], Set[str], List[str]]:
    """
    Read TSV into list of dict rows.
    Collect:
      - uniprots: primary_accession values
      - iprs: all InterPro IDs from 'interpro_ids' (split by ';')
    Return (rows, iprs, uniprots, header_fields).
    """
    rows: List[Dict[str,str]] = []
    iprs: Set[str] = set()
    uniprots: Set[str] = set()

    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        in_header = next(reader)

    total_rows = count_data_rows(tsv_path)
    if TQDM:
        pbar = _tqdm(total=total_rows, desc="[2/4] TSV rows", unit="row")
    else:
        _info("[2/4] Reading TSV rows…")

    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows.append(row)
            acc = (row.get("primary_accession") or "").strip()
            if acc:
                uniprots.add(acc)
            ipr_field = (row.get("interpro_ids") or "").strip()
            if ipr_field:
                for tok in ipr_field.split(";"):
                    ipr = tok.strip()
                    if ipr:
                        iprs.add(ipr)
            if TQDM:
                pbar.update(1)

    if TQDM:
        pbar.close()

    return rows, iprs, uniprots, in_header

def parse_protein2ipr_gz_to_map(gz_path: Path, needed_iprs: Set[str], needed_uniprots: Set[str], pigz_threads: int = 12) -> Dict[str, str]:
    """
    Stream protein2ipr.dat.gz as TSV with columns:
      0: UniProt accession   1: IPR ID   2: IPR description   (3+: ignored)
    Only consider rows where UniProt ∈ needed_uniprots AND IPR ∈ needed_iprs.
    Build map IPR -> first non-empty description.
    """
    ipr_to_desc: Dict[str, str] = {}
    first_line_seen = False
    lines_seen = 0

    if TQDM:
        pbar = _tqdm(total=None, desc="[3/4] protein2ipr.gz (stream)", unit="line", mininterval=0.2, smoothing=0.1)
        pbar.set_postfix_str("warming up…"); pbar.refresh()
    else:
        _info("[3/4] Parsing protein2ipr.dat.gz…")

    with open_text_auto(gz_path, pigz_threads=pigz_threads) as f:
        for raw in f:
            lines_seen += 1
            if TEST_LIMIT_GZ_LINES and lines_seen >= TEST_LIMIT_GZ_LINES:
                if TQDM:
                    pbar.set_postfix_str(f"line limit reached: {lines_seen}"); pbar.refresh()
                break

            if TQDM and not first_line_seen:
                first_line_seen = True
                pbar.set_postfix_str(""); pbar.refresh()

            if not raw or raw.startswith("#"):
                if TQDM: pbar.update(1)
                continue

            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 3:
                if TQDM: pbar.update(1)
                continue

            uni  = parts[0].strip()
            ipr  = parts[1].strip()
            desc = parts[2].strip()

            if not uni or not ipr:
                if TQDM: pbar.update(1)
                continue

            if needed_uniprots and uni not in needed_uniprots:
                if TQDM: pbar.update(1)
                continue
            if needed_iprs and ipr not in needed_iprs:
                if TQDM: pbar.update(1)
                continue

            if desc == "-" or desc is None:
                desc = ""

            prev = ipr_to_desc.get(ipr)
            if prev is None:
                ipr_to_desc[ipr] = desc
            elif not prev and desc:
                ipr_to_desc[ipr] = desc

            if TQDM: pbar.update(1)

    if TQDM:
        pbar.close()

    return ipr_to_desc

def format_interpro_with_ids_upper(ipr_ids_field: str, ipr_map: Dict[str,str]) -> str:
    """
    Given a semicolon-delimited list of raw IPR IDs (e.g., "IPR050125;IPR000276"),
    return "(IPR050125) DESCRIPTION; (IPR000276) DESCRIPTION; ..." (uppercase).
    If a description is missing, include "(IPRxxxxx);" with no desc.
    """
    if not ipr_ids_field:
        return ""
    out = []
    for raw in ipr_ids_field.split(";"):
        ipr = raw.strip().upper()
        if not ipr:
            continue
        desc = (ipr_map.get(ipr, "") or "").strip().upper()
        out.append(f"({ipr}) {desc};" if desc else f"({ipr});")
    return " ".join(out)

def main():
    in_path  = Path(INPUT_TSV)
    gz_path  = Path(PROTEIN2IPR_GZ)
    out_path = Path(OUTPUT_TSV)

    if not in_path.exists():
        _info(f"[ERROR] Input TSV not found: {in_path}")
        sys.exit(1)
    if not gz_path.exists():
        _info(f"[ERROR] protein2ipr.gz not found: {gz_path}")
        sys.exit(1)

    _info("[1/4] Counting TSV rows…")
    total_rows = count_data_rows(in_path)
    _info(f"[INFO] TSV data rows: {total_rows:,}")

    rows, needed_iprs, needed_uniprots, in_header = parse_tsv_collect(in_path)
    _info(f"[INFO] Unique InterPro IDs in TSV: {len(needed_iprs):,}")
    _info(f"[INFO] Unique UniProt accessions in TSV: {len(needed_uniprots):,}")

    # Normalize keys to uppercase so lookups work regardless of case
    needed_iprs_upper = {ipr.upper() for ipr in needed_iprs}
    needed_uniprots_upper = {u.upper() for u in needed_uniprots}

    ipr_map_raw = parse_protein2ipr_gz_to_map(gz_path, needed_iprs=needed_iprs_upper, needed_uniprots=needed_uniprots_upper, pigz_threads=12)
    # Normalize map keys to uppercase for consistent lookup
    ipr_map = {k.upper(): (v or "").upper() for k, v in ipr_map_raw.items()}
    _info(f"[INFO] Mapped IPR descriptions: {len(ipr_map):,}")

    out_fields = in_header  # keep original columns; we replace the contents of interpro_ids

    if TQDM:
        pbar = _tqdm(total=len(rows), desc="[4/4] Writing", unit="row")
    else:
        _info("[4/4] Writing output TSV…")

    with out_path.open("w", encoding="utf-8", newline="") as f_out:
        w = csv.DictWriter(f_out, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        for row in rows:
            # Uppercase protein and gene
            if "protein" in row and row["protein"] is not None:
                row["protein"] = row["protein"].strip().upper()
            if "gene" in row and row["gene"] is not None:
                row["gene"] = row["gene"].strip().upper()

            # Replace interpro_ids with uppercase formatted "(IPR) DESCRIPTION;"
            ipr_ids = (row.get("interpro_ids") or "").strip()
            row["interpro_ids"] = format_interpro_with_ids_upper(ipr_ids, ipr_map)
            w.writerow(row)
            if TQDM:
                pbar.update(1)

    if TQDM:
        pbar.close()

    _info(f"[DONE] Wrote: {out_path}")

if __name__ == "__main__":
    main()
