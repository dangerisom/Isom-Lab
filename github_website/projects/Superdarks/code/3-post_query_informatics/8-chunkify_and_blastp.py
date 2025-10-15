#!/usr/bin/env python3
"""
No-CLI chunk + BLASTP submitter (LSF bsub).

- Edit the CONFIG block below, then run:
      python 10-chunkify_and_blastp_nr.py
- Splits one FASTA into N chunks (whole-record).
- Submits one BLASTP job per chunk using LSF `bsub`.
- The BLASTP runner block below is the same as in script 8; we keep it
  unchanged for consistency. For cluster jobs, we build the equivalent
  blastp command and pass it directly to `bsub` (no command-line interface
  to this Python script).

Outputs (under OUT_DIRECTORY / <fasta_stem>):
  <OUT_DIRECTORY>/<stem>/chunks/*.fa
  <OUT_DIRECTORY>/<stem>/blastp_out/*.blastp.tsv
  <OUT_DIRECTORY>/<stem>/logs_blastp/*.%J.{out,err}
"""

from __future__ import annotations
import math
import shlex
import subprocess
import shutil
from pathlib import Path
from shutil import which

# =========================== CONFIG ===========================

# SINGLE input FASTA to split and BLAST:
FASTA_FILE  = Path("/nethome/dgi5/blastdb/ABE_reps.fa")

# BLAST database *prefix* (the path you used with makeblastdb -out ...), not .pin etc.
DB_NAME     = Path("/nethome/dgi5/blastdb/ABE_reps")

# ---------------------------------------------------------------

# # SINGLE input FASTA to split and BLAST:
# FASTA_FILE  = Path("/nethome/dgi5/blastdb/trim_blastp_sequences__Homo-sapiens.fa")

# # BLAST database *prefix* (the path you used with makeblastdb -out ...), not .pin etc.
# DB_NAME     = Path("/nethome/dgi5/blastdb/trim_blastp_sequences__Homo-sapiens")

# ---------------------------------------------------------------

# Where to place all outputs. We will create OUT_DIRECTORY / FASTA_FILE.stem / {chunks, blastp_out, logs_blastp}
OUT_DIRECTORY = Path("/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/2023.11.12.7tmps-analysis_final")

# Number of chunks to split FASTA into:
NUM_CHUNKS  = 256

# Optional LSF queue; set "" or None to use site default
QUEUE       = "isomlab"

# Per-job threads for BLASTP (passed to -num_threads). Tune for your cluster policy.
THREADS_PER_JOB = 1

# Path to YOUR blastp binary (absolute path recommended)
BLASTP_BIN  = "/nethome/dgi5/anaconda3/envs/latest_blast/bin/blastp"

# BLASTP parameters (kept aligned with script 8)
BLAST_MATRIX        = "BLOSUM62"   # use "BLOSUM45" to trigger special gap penalties
BLAST_EVALUE        = 1
BLAST_GAP_OPEN      = 11
BLAST_GAP_EXTEND    = 1
BLAST_SEG           = "yes"
BLAST_XDROP         = 15
BLAST_XDROP_FINAL   = 15
BLAST_MAX_TARGETS   = 1000
# =============================================================


# ------------------------ BLASTP runner (exactly as in script 8) ------------------------
def _resolve_exe(exe: str) -> str:
    """
    Return absolute path to the executable or raise if not found.
    Accepts absolute/relative path or a name resolvable via PATH.
    """
    p = Path(exe)
    if p.is_file():
        return str(p)
    w = which(exe)
    if w:
        return w
    raise FileNotFoundError(f"blastp binary not found: {exe}")

def run_blastp(query_file, db_name, out_file,
               matrix="BLOSUM62", e_value=1, gap_open=11, gap_extend=1,
               seg="yes", xdrop_gap=15, xdrop_gap_final=15, max_target_seqs=1000000,
               blastp_bin="blastp", num_threads=None):
    """
    Execute BLASTP with the given parameters. Saves results to 'out_file'.
    """
    exe = _resolve_exe(blastp_bin)

    if matrix == "BLOSUM45":
        cmd = [
            exe,
            "-query", query_file,
            "-db", db_name,
            "-out", out_file,
            "-matrix", matrix,
            "-evalue", str(e_value),
            "-gapopen", str(14),
            "-gapextend", str(2),
            "-seg", str(seg),
            "-xdrop_gap", str(xdrop_gap),
            "-xdrop_gap_final", str(xdrop_gap_final),
            "-max_target_seqs", str(max_target_seqs),
            "-outfmt", "7 qseqid pident length evalue bitscore stitle",
        ]
    else:
        cmd = [
            exe,
            "-query", query_file,
            "-db", db_name,
            "-out", out_file,
            "-matrix", matrix,
            "-evalue", str(e_value),
            "-gapopen", str(gap_open),
            "-gapextend", str(gap_extend),
            "-seg", str(seg),
            "-xdrop_gap", str(xdrop_gap),
            "-xdrop_gap_final", str(xdrop_gap_final),
            "-max_target_seqs", str(max_target_seqs),
            "-outfmt", "7 qseqid pident length evalue bitscore stitle",
        ]

    # Honor threads if provided
    if num_threads:
        cmd += ["-num_threads", str(num_threads)]

    try:
        print(f"Running BLASTp: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        print(f"Results saved to {out_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running BLASTp for {query_file}: {e}")
# --------------------------------------------------------------------------------------


# --------- Build-only version of the BLASTP command (for bsub payload) ---------
def build_blastp_cmd(query_file: str, db_name: str, out_file: str,
                     matrix: str = BLAST_MATRIX, e_value: float = BLAST_EVALUE,
                     gap_open: int = BLAST_GAP_OPEN, gap_extend: int = BLAST_GAP_EXTEND,
                     seg: str = BLAST_SEG, xdrop_gap: int = BLAST_XDROP,
                     xdrop_gap_final: int = BLAST_XDROP_FINAL,
                     max_target_seqs: int = BLAST_MAX_TARGETS,
                     blastp_bin: str = BLASTP_BIN, num_threads: int | None = None) -> list[str]:
    """
    Create the exact same BLASTP command line that run_blastp() would execute,
    but do NOT run it. We use this to hand to LSF `bsub`.
    """
    exe = _resolve_exe(blastp_bin)

    if matrix == "BLOSUM45":
        cmd = [
            exe,
            "-query", query_file,
            "-db", db_name,
            "-out", out_file,
            "-matrix", matrix,
            "-evalue", str(e_value),
            "-gapopen", str(14),
            "-gapextend", str(2),
            "-seg", str(seg),
            "-xdrop_gap", str(xdrop_gap),
            "-xdrop_gap_final", str(xdrop_gap_final),
            "-max_target_seqs", str(max_target_seqs),
            "-outfmt", "7 qseqid pident length evalue bitscore stitle",
        ]
    else:
        cmd = [
            exe,
            "-query", query_file,
            "-db", db_name,
            "-out", out_file,
            "-matrix", matrix,
            "-evalue", str(e_value),
            "-gapopen", str(gap_open),
            "-gapextend", str(gap_extend),
            "-seg", str(seg),
            "-xdrop_gap", str(xdrop_gap),
            "-xdrop_gap_final", str(xdrop_gap_final),
            "-max_target_seqs", str(max_target_seqs),
            "-outfmt", "7 qseqid pident length evalue bitscore stitle",
        ]

    if num_threads:
        cmd += ["-num_threads", str(num_threads)]
    return cmd


# ------------------------ FASTA chunking ------------------------
def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)


def count_fasta_sequences(fasta_path: Path) -> int:
    n = 0
    with fasta_path.open('r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.startswith('>'):
                n += 1
    return n


def chunk_fasta(in_fasta: Path, out_dir: Path, num_chunks: int) -> list[Path]:
    """
    Split FASTA into `num_chunks` files, preserving whole records.
    Returns list of chunk file paths in order.
    """
    if num_chunks < 1:
        raise ValueError("NUM_CHUNKS must be >= 1")
    ensure_dirs(out_dir)

    total = count_fasta_sequences(in_fasta)
    if total == 0:
        raise ValueError(f"No sequences found in {in_fasta}")

    per = math.ceil(total / num_chunks)

    stem = in_fasta.stem
    created: list[Path] = []
    writers = []
    try:
        # Open first chunk writer
        chunk_idx = 1
        chunk_path = out_dir / f"{stem}.chunk_{chunk_idx:05d}.fa"
        writers.append(chunk_path.open('w', encoding='utf-8'))
        created.append(chunk_path)
        written_in_current = 0

        with in_fasta.open('r', encoding='utf-8', errors='ignore') as fh:
            header = None
            seq = []
            for raw in fh:
                if raw.startswith('>'):
                    if header is not None:
                        if written_in_current >= per and chunk_idx < num_chunks:
                            writers[-1].close()
                            chunk_idx += 1
                            chunk_path = out_dir / f"{stem}.chunk_{chunk_idx:05d}.fa"
                            writers.append(chunk_path.open('w', encoding='utf-8'))
                            created.append(chunk_path)
                            written_in_current = 0
                        w = writers[-1]
                        w.write(header + '\n')
                        for s in seq:
                            w.write(s + '\n')
                        written_in_current += 1
                    header = raw.rstrip('\n')
                    seq = []
                else:
                    seq.append(raw.rstrip('\n'))
            # flush last
            if header is not None:
                if written_in_current >= per and chunk_idx < num_chunks:
                    writers[-1].close()
                    chunk_idx += 1
                    chunk_path = out_dir / f"{stem}.chunk_{chunk_idx:05d}.fa"
                    writers.append(chunk_path.open('w', encoding='utf-8'))
                    created.append(chunk_path)
                    written_in_current = 0
                w = writers[-1]
                w.write(header + '\n')
                for s in seq:
                    w.write(s + '\n')
                written_in_current += 1
    finally:
        for w in writers:
            try:
                w.close()
            except Exception:
                pass

    print(f"Chunked {total} sequences into {len(created)} files in {out_dir}")
    return created


# ------------------------ LSF job submission ------------------------
def bsub_submit(cmd: list[str], job_name: str, stdout_path: Path, stderr_path: Path, queue: str | None = None) -> None:
    """
    Submit a command to LSF via bsub. `cmd` is the actual program + args to run.
    """
    bsub = which("bsub")
    if not bsub:
        raise FileNotFoundError("bsub not found on PATH. Are you on your cluster login node?")

    bsub_cmd = [bsub, "-J", job_name, "-oo", str(stdout_path), "-eo", str(stderr_path)]
    if queue:
        bsub_cmd += ["-q", queue]

    # Robust quoting for the payload
    payload = " ".join(shlex.quote(x) for x in cmd)
    full = bsub_cmd + [payload]
    print("Submitting:", " ".join(full))
    subprocess.run(full, check=True)


def submit_blastp_jobs_for_chunks(chunks: list[Path],
                                  db_name: Path,
                                  out_dir: Path,
                                  logs_dir: Path,
                                  matrix: str = BLAST_MATRIX,
                                  e_value: float = BLAST_EVALUE,
                                  gap_open: int = BLAST_GAP_OPEN,
                                  gap_extend: int = BLAST_GAP_EXTEND,
                                  seg: str = BLAST_SEG,
                                  xdrop_gap: int = BLAST_XDROP,
                                  xdrop_gap_final: int = BLAST_XDROP_FINAL,
                                  max_target_seqs: int = BLAST_MAX_TARGETS,
                                  blastp_bin: str = BLASTP_BIN,
                                  num_threads: int = THREADS_PER_JOB,
                                  queue: str | None = QUEUE):
    """
    For each chunk FASTA, submit one BLASTP job. Output TSVs to out_dir, logs to logs_dir.
    Builds the exact BLAST command string (same logic as script 8) and hands it to bsub.
    """
    ensure_dirs(out_dir, logs_dir)

    for chunk in chunks:
        stem = chunk.stem  # e.g., myinput.chunk_00001
        out_tsv = out_dir / f"{stem}.blastp.tsv"
        job = f"blastp_{stem}"
        out_log = logs_dir / f"{stem}.%J.out"
        err_log = logs_dir / f"{stem}.%J.err"

        cmd = build_blastp_cmd(
            query_file=str(chunk),
            db_name=str(db_name),
            out_file=str(out_tsv),
            matrix=matrix,
            e_value=e_value,
            gap_open=gap_open,
            gap_extend=gap_extend,
            seg=seg,
            xdrop_gap=xdrop_gap,
            xdrop_gap_final=xdrop_gap_final,
            max_target_seqs=max_target_seqs,
            blastp_bin=blastp_bin,
            num_threads=num_threads,
        )
        bsub_submit(cmd, job, out_log, err_log, queue=queue)


# ------------------------ main (no CLI) ------------------------
def main():
    fasta = FASTA_FILE
    db = DB_NAME

    if not fasta.is_file():
        raise SystemExit(f"FASTA file not found: {fasta}")
    # Sanity note (not fatal): check presence of typical BLAST DB files.
    if not Path(str(db) + '.pin').exists() and not Path(str(db) + '.phr').exists():
        print(f"WARNING: I didn't find typical BLAST DB files next to {db}.")
        print("         Ensure DB_NAME is the prefix you used for makeblastdb -out.")

    # Base output directory derived from FASTA stem
    base_out = OUT_DIRECTORY / fasta.stem

    # --- Replace output directory if it already exists ---
    if base_out.exists():
        print(f"[INFO] Output directory exists; replacing: {base_out}")
        if base_out.is_dir():
            shutil.rmtree(base_out)
        else:
            base_out.unlink(missing_ok=True)

    chunks_dir = base_out / "chunks"
    out_dir    = base_out / "blastp_out"
    logs_dir   = base_out / "logs_blastp"
    ensure_dirs(chunks_dir, out_dir, logs_dir)

    # 1) Chunk the FASTA
    chunks = chunk_fasta(fasta, chunks_dir, NUM_CHUNKS)

    # 2) Submit one job per chunk
    submit_blastp_jobs_for_chunks(
        chunks=chunks,
        db_name=db,
        out_dir=out_dir,
        logs_dir=logs_dir,
        matrix=BLAST_MATRIX,
        e_value=BLAST_EVALUE,
        gap_open=BLAST_GAP_OPEN,
        gap_extend=BLAST_GAP_EXTEND,
        seg=BLAST_SEG,
        xdrop_gap=BLAST_XDROP,
        xdrop_gap_final=BLAST_XDROP_FINAL,
        max_target_seqs=BLAST_MAX_TARGETS,
        blastp_bin=BLASTP_BIN,
        num_threads=THREADS_PER_JOB,
        queue=(QUEUE or None),
    )

    print("\nSubmitted all jobs. Monitor with `bjobs` and check:", logs_dir)


if __name__ == "__main__":
    main()
