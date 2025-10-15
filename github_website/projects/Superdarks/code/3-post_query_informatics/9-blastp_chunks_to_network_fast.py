#!/usr/bin/env python3
"""
Unified network builder (FAST or TURBO/LSF) from BLASTP chunk outputs.

Outputs (compatible with scripts 12 / 13):
  network_out/nodes.tsv.gz            (uniprot_id, superkingdom, gene, organism, rank, degree_unique)
  network_out/edges_unique.tsv.gz     (src_uid, dst_uid, max_pident)
  network_out/edges_all.tsv.gz        (optional; src_uid, dst_uid, pident, evalue, length, bitscore)
  network_out/network_summary.txt

Assumes BLAST outfmt 7 (no headers inside file body):
  qseqid  pident  length  evalue  bitscore  stitle
and FASTA headers structured as:
  superkingdom|uniprot_id|gene_or_symbol|organism|rank
"""

from __future__ import annotations
import csv
import gzip
import math
import os
import re
import shlex
import sqlite3
import subprocess
import sys

def _resolve_python_bin(p: str) -> str:
    """Return absolute path to python binary if possible."""
    try:
        import shutil, os as _os
        if p and p.startswith("/") and Path(p).exists():
            return p
        found = shutil.which(p or "python")
        return found or (p or "python")
    except Exception:
        return p or "python"

import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from shutil import which
from typing import Dict, Iterable, List, Optional, Tuple

# ============================== CONFIG ===============================
# Where your per-chunk BLASTP TSVs live (files like *.blastp.tsv)
# BLASTP_DIR = Path("/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/2023.11.12.7tmps-analysis_final/2023.11.12.trim_blastp_sequences_final_nr_id0p50_c0p50/blastp_out")
# BLASTP_GLOB = "*.blastp.tsv"

# ---------------------------------------------------------------------

# BLASTP_DIR = Path("/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/2023.11.12.7tmps-analysis_final/ABE_reps/blastp_out")
# BLASTP_GLOB = "*.blastp.tsv"

# ---------------------------------------------------------------------

BLASTP_DIR = Path("/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/2023.11.12.7tmps-analysis_final/trim_blastp_sequences_human_gpcrs/blastp_out")
BLASTP_GLOB = "*.blastp.tsv"

# Final outputs root (the script will create <OUT_ROOT>/network_out/*)
OUT_ROOT = BLASTP_DIR
NETWORK_OUT_DIR = OUT_ROOT / "network_out_fast"

# Keep a huge per-hit file? (src dst pident evalue length bitscore)
WRITE_EDGES_ALL = True

# Optional early filter: drop rows with pident < PIDENT_MIN while building.
# (Leave at 0.0 to keep everything and let 12/13 do their own thresholding.)
PIDENT_MIN = 0.0

# ----------------- MODE picker -----------------
# "fast"  = single-node, bulk-SQL optimized
# "turbo" = launch LSF jobs to build per-worker partials, then reduce
MODE = "turbo"

# TURBO (LSF) settings
N_WORKERS = 64                  # array width / number of separate jobs
QUEUE = "isomlab"               # set "" to use default site queue
PYTHON_BIN = sys.executable     # python path used by workers
WAIT_FOR_WORKERS = True         # coordinator polls for completion
POLL_SECS = 30                  # how often to poll for completion

# Batched submission to avoid stampeding the scheduler
BATCH_SUBMIT_SIZE = int(os.environ.get("BATCH_SUBMIT_SIZE", "8"))

# Allow env overrides for convenience
N_WORKERS = int(os.environ.get("N_WORKERS", str(N_WORKERS)))
QUEUE = os.environ.get("QUEUE", QUEUE)
PYTHON_BIN = os.environ.get("PYTHON_BIN", PYTHON_BIN)
WAIT_FOR_WORKERS = os.environ.get("WAIT_FOR_WORKERS", str(WAIT_FOR_WORKERS)).lower() in ("1","true","yes","y")
POLL_SECS = int(os.environ.get("POLL_SECS", str(POLL_SECS)))
# ====================================================================

# Optional LSF knobs (override via env)
LSF_PROJECT  = os.environ.get("LSF_PROJECT", "")        # e.g., "isomlab_proj"
LSF_RES      = os.environ.get("LSF_RES", "")            # e.g., 'rusage[mem=4000]'
LSF_RUNTIME  = os.environ.get("LSF_RUNTIME", "04:00")   # walltime HH:MM
LSF_CORES    = int(os.environ.get("LSF_CORES", "1"))    # slots per worker
LSF_MEMORYMB = int(os.environ.get("LSF_MEMORYMB", "4000"))
LSF_STDOUT_EXT = ".%J.out"
LSF_STDERR_EXT = ".%J.err"


# ----------- header parsing helpers -----------
HEADER_SPLIT_RE = re.compile(r"\|")
WS_RE = re.compile(r"\s+", flags=re.UNICODE)

@dataclass(frozen=True)
class NodeAttrs:
    superkingdom: str
    uid: str
    gene: str
    organism: str
    rank: str

def _sanitize_token(s: str) -> str:
    if s is None:
        return ""
    s = s.strip()
    s = s.replace("\u00A0", " ").replace("\u2007", " ").replace("\u202F", " ")
    s = WS_RE.sub("_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s

def parse_header_to_attrs(header: str) -> Optional[NodeAttrs]:
    """
    Eukaryota|A0A2U4ZXY9|COI|Tenthredo_sp._BBHYL871-10|2
    -> NodeAttrs
    """
    if not header:
        return None
    parts = HEADER_SPLIT_RE.split(header.strip())
    if len(parts) < 5:
        return None
    sk, uid, gene, org, rank = parts[0], parts[1], parts[2], parts[3], parts[4]
    return NodeAttrs(_sanitize_token(sk),
                     _sanitize_token(uid),
                     _sanitize_token(gene) or "UC",
                     _sanitize_token(org),
                     _sanitize_token(rank))


# ----------- I/O / SQLite utilities -----------
def ensure_clean_dir(p: Path) -> None:
    if p.exists():
        # rm -rf p
        for child in sorted(p.rglob("*"), reverse=True):
            try:
                child.unlink()
            except IsADirectoryError:
                try: child.rmdir()
                except Exception: pass
        try: p.rmdir()
        except Exception: pass
    p.mkdir(parents=True, exist_ok=True)

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def sqlite_connect(db_path: Path) -> sqlite3.Connection:
    con = sqlite3.connect(str(db_path))
    con.execute("PRAGMA journal_mode=WAL;")
    con.execute("PRAGMA synchronous=OFF;")
    con.execute("PRAGMA temp_store=MEMORY;")
    con.execute("PRAGMA mmap_size=30000000000;")
    con.execute("PRAGMA cache_size=-100000;")  # ~100MB
    return con

def iter_blast_rows(tsv_path: Path):
    """
    Yields: (qseqid, pident, length, evalue, bitscore, stitle)
    """
    with tsv_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or (row and row[0].startswith("#")):
                continue
            if len(row) < 6:
                continue
            qseqid = row[0].strip()
            try:
                pident = float(row[1])
                length = int(row[2])
                evalue = float(row[3]) if row[3] not in ("", "NA") else math.inf
                bits   = float(row[4])
            except Exception:
                continue
            stitle = row[5].strip()
            yield (qseqid, pident, length, evalue, bits, stitle)


# =============== FAST path (single node, bulk SQL) ===============
def fast_build(files: List[Path], out_dir: Path) -> None:
    tmp_dir = out_dir / "tmp_fast"
    ensure_clean_dir(tmp_dir)
    ensure_dir(out_dir)

    # optional big "all edges" file
    eall = None
    if WRITE_EDGES_ALL:
        eall = gzip.open(out_dir / "edges_all.tsv.gz", "wt", encoding="utf-8")
        eall.write("src_uid\tdst_uid\tpident\tevalue\tlength\tbitscore\n")

    # Collect node attributes (first seen wins)
    nodes: Dict[str, NodeAttrs] = {}

    # Stage raw edges in SQLite, then aggregate to unique pairs in one SQL
    db = tmp_dir / "stage.sqlite3"
    con = sqlite_connect(db)
    cur = con.cursor()

    cur.execute("""
        CREATE TABLE stg_edges (
            a TEXT NOT NULL,
            b TEXT NOT NULL,
            pident REAL NOT NULL,
            evalue REAL NOT NULL,
            bitscore REAL NOT NULL,
            length INTEGER NOT NULL
        );
    """)
    # covering index for group by
    cur.execute("CREATE INDEX idx_stg_ab ON stg_edges(a,b);")
    con.commit()

    BATCH = 100000
    buf: List[Tuple[str,str,float,float,float,int]] = []
    total_rows = 0
    kept_rows = 0
    bad_lines = 0

    for f in files:
        for qseqid, pident, length, evalue, bits, stitle in iter_blast_rows(f):
            total_rows += 1
            if pident < PIDENT_MIN:
                continue
            qa = parse_header_to_attrs(qseqid)
            sa = parse_header_to_attrs(stitle)
            if qa is None or sa is None or not qa.uid or not sa.uid:
                bad_lines += 1
                continue

            if qa.uid not in nodes: nodes[qa.uid] = qa
            if sa.uid not in nodes: nodes[sa.uid] = sa

            a, b = (qa.uid, sa.uid) if qa.uid <= sa.uid else (sa.uid, qa.uid)
            if a == b:
                continue

            if eall:
                eall.write(f"{qa.uid}\t{sa.uid}\t{pident:.3f}\t{(evalue if math.isfinite(evalue) else 'inf')}\t{length}\t{bits:.1f}\n")

            buf.append((a, b, pident, (evalue if math.isfinite(evalue) else 1e308), bits, length))
            kept_rows += 1

            if len(buf) >= BATCH:
                cur.executemany("INSERT INTO stg_edges(a,b,pident,evalue,bitscore,length) VALUES (?,?,?,?,?,?)", buf)
                con.commit()
                buf.clear()

    if buf:
        cur.executemany("INSERT INTO stg_edges(a,b,pident,evalue,bitscore,length) VALUES (?,?,?,?,?,?)", buf)
        con.commit()
        buf.clear()

    if eall:
        eall.close()

    # Aggregate to unique undirected edges with best stats
    cur.execute("""
        CREATE TABLE edges_agg (
            a TEXT NOT NULL,
            b TEXT NOT NULL,
            pident_max REAL NOT NULL,
            evalue_min REAL NOT NULL,
            bitscore_max REAL NOT NULL,
            length_max INTEGER NOT NULL,
            PRIMARY KEY (a,b)
        ) WITHOUT ROWID;
    """)
    con.commit()

    cur.execute("""
        INSERT INTO edges_agg(a,b,pident_max,evalue_min,bitscore_max,length_max)
        SELECT a, b,
               MAX(pident)          AS pident_max,
               MIN(evalue)          AS evalue_min,
               MAX(bitscore)        AS bitscore_max,
               MAX(length)          AS length_max
        FROM stg_edges
        GROUP BY a,b;
    """)
    con.commit()

    # Degrees
    degree: Counter[str] = Counter()
    for uid, cnt in cur.execute("SELECT a, COUNT(*) FROM edges_agg GROUP BY a;"):
        degree[uid] += int(cnt)
    for uid, cnt in cur.execute("SELECT b, COUNT(*) FROM edges_agg GROUP BY b;"):
        degree[uid] += int(cnt)

    # Write edges_unique.tsv.gz (6 columns: src_uid dst_uid pident evalue length bitscore)
    with gzip.open(out_dir / "edges_unique.tsv.gz", "wt", encoding="utf-8") as eu:
        eu.write("src_uid\tdst_uid\tpident\tevalue\tlength\tbitscore\n")
        for a, b, pid, emin, lmax, bmax in cur.execute("SELECT a,b,pident_max,evalue_min,length_max,bitscore_max FROM edges_agg ORDER BY a,b;"):
            eu.write(f"{a}\t{b}\t{pid:.3f}\t{(emin if emin != 1e308 else float('inf'))}\t{lmax}\t{bmax:.1f}\n")

    # Write nodes.tsv.gz
    with gzip.open(out_dir / "nodes.tsv.gz", "wt", encoding="utf-8") as nf:
        nf.write("uniprot_id\tsuperkingdom\tgene\torganism\trank\tdegree_unique\n")
        for uid, attrs in nodes.items():
            nf.write(f"{uid}\t{attrs.superkingdom}\t{attrs.gene}\t{attrs.organism}\t{attrs.rank}\t{degree.get(uid,0)}\n")

    # Summary
    with (out_dir / "network_summary.txt").open("w", encoding="utf-8") as sf:
        n_edges_unique = cur.execute("SELECT COUNT(*) FROM edges_agg").fetchone()[0]
        sf.write("=== Network Summary (FAST) ===\n")
        sf.write(f"BLAST dir           : {BLASTP_DIR}\n")
        sf.write(f"Files processed     : {len(files)}\n")
        sf.write(f"Total rows seen     : {total_rows:,}\n")
        sf.write(f"Kept rows (>=pid)   : {kept_rows:,}  (pid>={PIDENT_MIN})\n")
        sf.write(f"Bad/skip lines      : {bad_lines:,}\n")
        sf.write(f"Unique edges        : {n_edges_unique:,}\n")
        sf.write(f"Nodes               : {len(nodes):,}\n")
        if WRITE_EDGES_ALL:
            sf.write(f"Edges-all           : {out_dir/'edges_all.tsv.gz'}\n")

    cur.close()
    con.close()


# =============== TURBO path (LSF parallel) ===============


def bsub_submit(cmd: list[str], job_name: str, stdout_path: Path, stderr_path: Path,
                queue: str | None = None, cwd: Path | None = None) -> None:
    """
    Submit a command to LSF via bsub using a short payload:
    - create worker_dir/run_worker.sh that sets env and calls python
    - submit: bsub [opts] -cwd worker_dir /bin/bash -lc ./run_worker.sh
    """
    bsub = which("bsub")
    if not bsub:
        raise FileNotFoundError("bsub not found on PATH. Are you on your cluster login node?")

    # Ensure directories & absolute paths
    cwd_abs = cwd.resolve() if cwd else Path.cwd().resolve()
    cwd_abs.mkdir(parents=True, exist_ok=True)
    stdout_abs = stdout_path.resolve()
    stderr_abs = stderr_path.resolve()
    stdout_abs.parent.mkdir(parents=True, exist_ok=True)
    stderr_abs.parent.mkdir(parents=True, exist_ok=True)

    # Build a small wrapper script to avoid very long payloads
    payload_str = " ".join(shlex.quote(x) for x in cmd)
    wrapper = cwd_abs / "run_worker.sh"
    script_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "echo \"[wrapper] host=$(hostname) user=$(whoami) pwd=$(pwd)\"",
        f"echo \"[wrapper] running: {payload_str}\"",
        payload_str,
    ]
    wrapper.write_text("\n".join(script_lines) + "\n", encoding="utf-8")
    wrapper.chmod(0o755)

    # Compose bsub command: short payload
    bsub_cmd = [bsub, "-J", job_name, "-oo", str(stdout_abs), "-eo", str(stderr_abs)]
    if queue:
        bsub_cmd += ["-q", queue]
    bsub_cmd += ["-cwd", str(cwd_abs), "/bin/bash", "-lc", "./run_worker.sh"]

    print("Submitting:", " ".join(bsub_cmd))
    subprocess.run(bsub_cmd, check=True)


def split_files_round_robin(files: List[Path], n: int) -> List[List[Path]]:
    buckets = [[] for _ in range(n)]
    for i, f in enumerate(sorted(files)):
        buckets[i % n].append(f)
    return buckets

def worker_process(manifest: Path, worker_dir: Path) -> None:
    """
    Each worker:
      - reads its file list from manifest
      - stages to sqlite stg_edges
      - aggregates to edges_agg
      - writes:
          worker_edges.tsv.gz   (a,b,pident_max,evalue_min,bitscore_max,length_max)
          worker_nodes.tsv.gz   (uid, superkingdom, gene, organism, rank)
      - touches done.ok
    """
    ensure_dir(worker_dir)
    db = worker_dir / "worker.sqlite3"
    con = sqlite_connect(db)
    cur = con.cursor()
    cur.execute("""
        CREATE TABLE stg_edges (
            a TEXT NOT NULL,
            b TEXT NOT NULL,
            pident REAL NOT NULL,
            evalue REAL NOT NULL,
            bitscore REAL NOT NULL,
            length INTEGER NOT NULL
        );
    """)
    cur.execute("CREATE INDEX idx_stg_ab ON stg_edges(a,b);")
    con.commit()

    nodes: Dict[str, NodeAttrs] = {}
    BATCH = 100000
    buf: List[Tuple[str,str,float,float,float,int]] = []

    files = [Path(x.strip()) for x in manifest.read_text(encoding="utf-8").splitlines() if x.strip()]
    kept_rows = 0
    for f in files:
        for qseqid, pident, length, evalue, bits, stitle in iter_blast_rows(f):
            if pident < PIDENT_MIN:
                continue
            qa = parse_header_to_attrs(qseqid)
            sa = parse_header_to_attrs(stitle)
            if qa is None or sa is None or not qa.uid or not sa.uid:
                continue
            if qa.uid not in nodes: nodes[qa.uid] = qa
            if sa.uid not in nodes: nodes[sa.uid] = sa
            a, b = (qa.uid, sa.uid) if qa.uid <= sa.uid else (sa.uid, qa.uid)
            if a == b:
                continue
            buf.append((a, b, pident, (evalue if math.isfinite(evalue) else 1e308), bits, length))
            kept_rows += 1
            if len(buf) >= BATCH:
                cur.executemany("INSERT INTO stg_edges(a,b,pident,evalue,bitscore,length) VALUES (?,?,?,?,?,?)", buf)
                con.commit()
                buf.clear()

    if buf:
        cur.executemany("INSERT INTO stg_edges(a,b,pident,evalue,bitscore,length) VALUES (?,?,?,?,?,?)", buf)
        con.commit()
        buf.clear()

    cur.execute("""
        CREATE TABLE edges_agg (
            a TEXT NOT NULL,
            b TEXT NOT NULL,
            pident_max REAL NOT NULL,
            evalue_min REAL NOT NULL,
            bitscore_max REAL NOT NULL,
            length_max INTEGER NOT NULL,
            PRIMARY KEY (a,b)
        ) WITHOUT ROWID;
    """)
    con.commit()

    cur.execute("""
        INSERT INTO edges_agg(a,b,pident_max,evalue_min,bitscore_max,length_max)
        SELECT a, b,
               MAX(pident),
               MIN(evalue),
               MAX(bitscore),
               MAX(length)
        FROM stg_edges
        GROUP BY a,b;
    """)
    con.commit()

    # dump worker products
    with gzip.open(worker_dir / "worker_edges.tsv.gz", "wt", encoding="utf-8") as ef:
        # a,b,pident_max,evalue_min,bitscore_max,length_max
        for a, b, pid, ev, bs, ln in cur.execute("SELECT a,b,pident_max,evalue_min,bitscore_max,length_max FROM edges_agg ORDER BY a,b;"):
            ef.write(f"{a}\t{b}\t{pid:.3f}\t{ev:.6g}\t{bs:.1f}\t{ln}\n")

    with gzip.open(worker_dir / "worker_nodes.tsv.gz", "wt", encoding="utf-8") as nf:
        nf.write("uid\tsuperkingdom\tgene\torganism\trank\n")
        for uid, attrs in nodes.items():
            nf.write(f"{uid}\t{attrs.superkingdom}\t{attrs.gene}\t{attrs.organism}\t{attrs.rank}\n")

    Path(worker_dir / "done.ok").write_text("ok", encoding="utf-8")

    cur.close()
    con.close()

def reducer_collect(workers_root: Path, out_dir: Path) -> None:
    """
    Reduce worker outputs:
      - union worker_edges.tsv.gz into stg_edges, aggregate -> edges_agg
      - union worker_nodes.tsv.gz first-seen wins
      - compute degrees, write final outputs
    """
    ensure_dir(out_dir)
    tmp = out_dir / "tmp_reduce"
    ensure_clean_dir(tmp)
    db = tmp / "reduce.sqlite3"
    con = sqlite_connect(db)
    cur = con.cursor()

    cur.execute("""
        CREATE TABLE stg_edges (
            a TEXT NOT NULL,
            b TEXT NOT NULL,
            pident_max REAL NOT NULL,
            evalue_min REAL NOT NULL,
            bitscore_max REAL NOT NULL,
            length_max INTEGER NOT NULL
        );
    """)
    cur.execute("CREATE INDEX idx_stg_ab ON stg_edges(a,b);")
    con.commit()

    # Load worker edge parts
    BATCH = 200000
    buf = []
    n_parts = 0
    for wdir in sorted(workers_root.glob("worker_*")):
        part = wdir / "worker_edges.tsv.gz"
        if not part.exists():
            continue
        n_parts += 1
        with gzip.open(part, "rt", encoding="utf-8") as fh:
            for line in fh:
                a, b, pid, ev, bs, ln = line.rstrip("\n").split("\t")
                buf.append((a, b, float(pid), float(ev), float(bs), int(ln)))
                if len(buf) >= BATCH:
                    cur.executemany("INSERT INTO stg_edges(a,b,pident_max,evalue_min,bitscore_max,length_max) VALUES (?,?,?,?,?,?)", buf)
                    con.commit()
                    buf.clear()
    if buf:
        cur.executemany("INSERT INTO stg_edges(a,b,pident_max,evalue_min,bitscore_max,length_max) VALUES (?,?,?,?,?,?)", buf)
        con.commit()
        buf.clear()

    # Final aggregate (NOTE: this aggregates across workers)
    cur.execute("""
        CREATE TABLE edges_agg (
            a TEXT NOT NULL,
            b TEXT NOT NULL,
            pident_max REAL NOT NULL,
            evalue_min REAL NOT NULL,
            bitscore_max REAL NOT NULL,
            length_max INTEGER NOT NULL,
            PRIMARY KEY (a,b)
        ) WITHOUT ROWID;
    """)
    con.commit()
    cur.execute("""
        INSERT INTO edges_agg(a,b,pident_max,evalue_min,bitscore_max,length_max)
        SELECT a, b,
               MAX(pident_max),
               MIN(evalue_min),
               MAX(bitscore_max),
               MAX(length_max)
        FROM stg_edges
        GROUP BY a,b;
    """)
    con.commit()

    # Degrees
    degree: Counter[str] = Counter()
    for uid, cnt in cur.execute("SELECT a, COUNT(*) FROM edges_agg GROUP BY a;"):
        degree[uid] += int(cnt)
    for uid, cnt in cur.execute("SELECT b, COUNT(*) FROM edges_agg GROUP BY b;"):
        degree[uid] += int(cnt)

    # Merge nodes (first seen wins)
    nodes: Dict[str, NodeAttrs] = {}
    for wdir in sorted(workers_root.glob("worker_*")):
        part = wdir / "worker_nodes.tsv.gz"
        if not part.exists():
            continue
        with gzip.open(part, "rt", encoding="utf-8") as fh:
            header = next(fh, None)
            for line in fh:
                uid, sk, gene, org, rank = line.rstrip("\n").split("\t")
                if uid not in nodes:
                    nodes[uid] = NodeAttrs(sk, uid, gene, org, rank)

    # Final files
    with gzip.open(out_dir / "edges_unique.tsv.gz", "wt", encoding="utf-8") as eu:
        eu.write("src_uid\tdst_uid\tpident\tevalue\tlength\tbitscore\n")
        for a, b, pid, emin, lmax, bmax in cur.execute("SELECT a,b,pident_max,evalue_min,length_max,bitscore_max FROM edges_agg ORDER BY a,b;"):
            eu.write(f"{a}\t{b}\t{pid:.3f}\t{(emin if emin != 1e308 else float('inf'))}\t{lmax}\t{bmax:.1f}\n")

    if WRITE_EDGES_ALL:
        # Not reproduced in TURBO to save I/O; emit a tiny stub to document choice
        with gzip.open(out_dir / "edges_all.tsv.gz", "wt", encoding="utf-8") as fh:
            fh.write("# omitted in TURBO reduce to save I/O; set WRITE_EDGES_ALL=False or use FAST mode to produce per-hit file\n")

    with gzip.open(out_dir / "nodes.tsv.gz", "wt", encoding="utf-8") as nf:
        nf.write("uniprot_id\tsuperkingdom\tgene\torganism\trank\tdegree_unique\n")
        for uid, attrs in nodes.items():
            nf.write(f"{uid}\t{attrs.superkingdom}\t{attrs.gene}\t{attrs.organism}\t{attrs.rank}\t{degree.get(uid,0)}\n")

    with (out_dir / "network_summary.txt").open("w", encoding="utf-8") as sf:
        n_edges_unique = cur.execute("SELECT COUNT(*) FROM edges_agg").fetchone()[0]
        sf.write("=== Network Summary (TURBO reduce) ===\n")
        sf.write(f"Workers parts loaded: {n_parts}\n")
        sf.write(f"Unique edges        : {n_edges_unique:,}\n")
        sf.write(f"Nodes               : {len(nodes):,}\n")

    cur.close()
    con.close()



def _bjobs_active(job_ids: List[str]) -> int:
    """Return count of active jobs among job_ids (RUN or PEND). If bjobs missing, return -1."""
    bjobs = which("bjobs")
    if not bjobs or not job_ids:
        return -1
    try:
        res = subprocess.run([bjobs, "-noheader"] + job_ids, capture_output=True, text=True, check=False)
        out = (res.stdout or "") + "\n" + (res.stderr or "")
        active = 0
        for line in out.splitlines():
            if " RUN " in line or " PEND " in line:
                active += 1
        return active
    except Exception:
        return -1

def _tail_file(path: Path, n: int = 50) -> str:
    try:
        if not path.exists():
            return "(no file)"
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
        return "\n".join(lines[-n:])
    except Exception as e:
        return f"(error reading {path}: {e})"



def _job_entries(jobs):
    """Yield (i, wdir, jid) for each job tuple that may be (i,wdir) or (i,wdir,jid)."""
    for t in jobs:
        if len(t) == 3:
            i, wdir, jid = t
        elif len(t) == 2:
            i, wdir = t
            jid = None
        else:
            # unexpected shape; try best-effort
            try:
                i = t[0]; wdir = t[1]; jid = t[2] if len(t) > 2 else None
            except Exception:
                continue
        yield i, wdir, jid



def _launch_local_worker(manifest: Path, wdir: Path) -> None:
    """Run a worker locally as a blocking subprocess (fallback when LSF won't start jobs)."""
    env = os.environ.copy()
    env["UNIFIED_MODE"] = "TURBO_WORKER"
    env["WORKER_MANIFEST"] = str(manifest)
    env["WORKER_DIR"] = str(wdir)
    env["PIDENT_MIN"] = str(PIDENT_MIN)
    wdir.mkdir(parents=True, exist_ok=True)
    (wdir / "local_fallback.txt").write_text("running locally", encoding="utf-8")
    cmd = [PYTHON_BIN, str(Path(__file__).resolve())]
    print("[fallback] launching local worker:", " ".join(cmd), "cwd=", str(wdir))
    res = subprocess.run(cmd, cwd=str(wdir), env=env)
    if res.returncode != 0:
        raise SystemExit(f"Local worker failed with code {res.returncode} for {wdir}")


def turbo_coordinator(files: List[Path], out_dir: Path) -> None:
    ensure_dir(out_dir)
    print(f"[paths] OUT_ROOT={OUT_ROOT.resolve()}")
    print(f"[paths] NETWORK_OUT_DIR={out_dir.resolve()}")
    workers_root = out_dir / "workers"
    manifests_dir = out_dir / "manifests"
    logs_dir = out_dir / "logs"
    ensure_clean_dir(workers_root)
    ensure_clean_dir(manifests_dir)
    ensure_dir(logs_dir)

    buckets = split_files_round_robin(files, N_WORKERS)

    # Write manifests and submit one job per worker (simple, robust)
    script_path = Path(__file__).resolve()
    jobs = []
    for i, flist in enumerate(buckets, start=1):
        if not flist:
            continue
        manifest = manifests_dir / f"worker_{i:04d}.list"
        manifest.write_text("\n".join(str(x) for x in flist) + "\n", encoding="utf-8")
        wdir = workers_root / f"worker_{i:04d}"
        ensure_dir(wdir)
        log_base = logs_dir / f"worker_{i:04d}"
        env_kv = {
            "UNIFIED_MODE": "TURBO_WORKER",
            "WORKER_MANIFEST": str(manifest),
            "WORKER_DIR": str(wdir),
            "PIDENT_MIN": str(PIDENT_MIN),
        }
        cmd = ["/usr/bin/env"] + [f"{k}={v}" for k,v in env_kv.items()] + [str(PYTHON_BIN), str(script_path)]
        job_name = f"netw_w{i:04d}"
        stdout_path = logs_dir / f"worker_{i:04d}.%J.out"
        stderr_path = logs_dir / f"worker_{i:04d}.%J.err"
        bsub_submit(cmd, job_name, stdout_path, stderr_path, queue=(QUEUE or None), cwd=wdir)
        jobs.append((i, wdir))

    if not WAIT_FOR_WORKERS:
        print("[info] Submitted workers and exiting (WAIT_FOR_WORKERS=False).")
        return

    # Wait for done.ok files, and detect if nothing is running anymore
    print(f"[wait] Waiting for {len(jobs)} workers to finish ...")
    remaining = {i for i,_,_ in _job_entries(jobs)}
    job_id_map = {i: jid for i,_,jid in _job_entries(jobs) if jid}
    # Immediate check in case some finished very fast
    pre_done = []
    for i, wdir, _jid in _job_entries(jobs):
        if i in remaining and (wdir / "done.ok").is_file():
            pre_done.append(i)
    for i in pre_done:
        remaining.discard(i)
        print(f"[wait] worker {i:04d} already done")
    idle_cycles = 0
    while remaining:
        print("[wait] remaining:", ", ".join(f"{i:04d}" for i in sorted(remaining)))
        time.sleep(POLL_SECS)
        done_now = []
        for i, wdir, _jid in _job_entries(jobs):
            if i in remaining and (wdir / "done.ok").is_file():
                done_now.append(i)
        for i in done_now:
            remaining.discard(i)
            print(f"[wait] worker {i:04d} done")
        active = _bjobs_active([job_id_map.get(i) for i in remaining if job_id_map.get(i)])
        if active == 0 and remaining:
            idle_cycles += 1
        else:
            idle_cycles = 0
        if idle_cycles >= 2 and remaining:
            print("[warn] No active LSF jobs; launching remaining workers locally...")
            for i, wdir, _jid in _job_entries(jobs):
                if i in remaining and not (wdir / "done.ok").is_file():
                    # Reconstruct manifest path
                    manifest = (wdir.parent.parent / "manifests" / f"worker_{i:04d}.list")
                    _launch_local_worker(manifest, wdir)
                    if (wdir / "done.ok").is_file():
                        print(f"[wait] worker {i:04d} done (local fallback)")
                        remaining.discard(i)
            # continue loop; if something still stuck, logs will print next cycles
    print("[wait] All workers done. Reducing...")

    reducer_collect(workers_root, out_dir)


# =============== main ===============
def main():
    # Remove previous output if it exists
    ensure_clean_dir(NETWORK_OUT_DIR)
    files = sorted(BLASTP_DIR.glob(BLASTP_GLOB))
    if not files:
        raise SystemExit(f"No files matched {BLASTP_GLOB} in {BLASTP_DIR}")

    # Worker branch? (env-driven, no CLI needed)
    if os.getenv("UNIFIED_MODE", "") == "TURBO_WORKER":
        manifest = Path(os.environ["WORKER_MANIFEST"])
        worker_dir = Path(os.environ["WORKER_DIR"])
        (worker_dir / "started.ok").write_text("started", encoding="utf-8")

        # Optional: allow coordinator to bump threshold via env
        try:
            pm = float(os.getenv("PIDENT_MIN", str(PIDENT_MIN)))
            globals()["PIDENT_MIN"] = pm
        except Exception:
            pass
        print(f"[worker] manifest={manifest} -> {worker_dir}  (pid>={PIDENT_MIN})")
        worker_process(manifest, worker_dir)
        print("[worker] done.")
        return

    # Coordinator branch
    print(f"[mode] {MODE.upper()}  files={len(files)}  pid>={PIDENT_MIN}  edges_all={'yes' if WRITE_EDGES_ALL else 'no'}")
    out_dir = NETWORK_OUT_DIR
    if MODE.lower() == "fast":
        fast_build(files, out_dir)
        print("[FAST] done.")
    elif MODE.lower() == "turbo":
        turbo_coordinator(files, out_dir)
        print("[TURBO] done.")
    else:
        raise SystemExit(f"Unknown MODE={MODE!r} (use 'fast' or 'turbo').")


if __name__ == "__main__":
    # Worker-first entrypoint: when launched by LSF, do NOT run main()
    if os.getenv("UNIFIED_MODE", "") == "TURBO_WORKER":
        manifest = Path(os.environ["WORKER_MANIFEST"])
        worker_dir = Path(os.environ["WORKER_DIR"])
        try:
            (worker_dir / "started.ok").write_text("started", encoding="utf-8")
        except Exception:
            pass
        try:
            worker_process(manifest, worker_dir)
        except Exception:
            import traceback
            (worker_dir / "worker_internal.err").write_text(traceback.format_exc(), encoding="utf-8")
            raise
        finally:
            # Ensure a done.ok is written if worker completed without raising inside turbo_worker
            if (worker_dir / "done.ok").exists() is False:
                # turbo_worker should create it; if not, create as last resort
                try:
                    (worker_dir / "done.ok").write_text("done", encoding="utf-8")
                except Exception:
                    pass
        sys.exit(0)
    # Coordinator mode
    main()
