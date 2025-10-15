
"""
fasta_filter_and_pair_pipeline.py  (pipe-2nd token with UniProt fallback)
------------------------------------------------------------------------
Header parsing:
  - If the first whitespace-delimited token after '>' contains '|', we take
    the **second** '|' token as the UniProt accession (user rule).
  - Otherwise (or if it looks empty), we fall back to a UniProt regex that
    matches common 6-char and 10-char accessions anywhere in the header:
      * 6-char: e.g., Q9H9K5 (OPQ… formats)
      * 10-char: e.g., A0A0D0SLZ3 (A0A… formats)

Pipeline:
  1) Filter input FASTA by a UniProt ID list.
  2) Write overlapping pair FASTAs in the order of that ID list.
"""

import os
import re
import gzip
import subprocess
from pathlib import Path
# ---- Optional PDB parser import guard (for tar->PDB step) ----
try:
    from pdbFile import PDBfile as _PDBfile
    _PDBFILE_AVAILABLE = True
except Exception:
    _PDBFILE_AVAILABLE = False
# -------------------------------------------------------------
import shutil
from typing import Iterable, Set, Dict, Tuple, List

# ---------------------------- CONFIG
PDB_DIR            = r"/Users/danisom/Desktop/local_informatics/0-inputs/hrh1_traceback"  # If set, prefer this over PDB_TARBALL (edit me) ----------------------------
IN_FASTA            = r"/Users/danisom/Desktop/local_informatics/0-sequences/2023.11.12.trim_blastp_sequences.fa"  # e.g., r"C:\data\huge_sequences.fasta.gz" or "/Volumes/Disk/seqs.fasta.gz"
OUT_FILTERED_FASTA  = r"/Users/danisom/Desktop/local_informatics/0-sequences/hrh1_traceback.fa"  # e.g., r"C:\data\filtered_sequences.fa"
WANTED_IDS_PATH     = r"/Users/danisom/Desktop/local_informatics/0-inputs/hrh1_traceback.txt"  # text file with UniProt IDs (one per line or comma-separated)
OUT_PAIRS_DIR       = r"/Users/danisom/Desktop/local_informatics/0-sequences/hrh1_traceback"  # e.g., r"C:\data\pairs_out"
PDB_TARBALL          = r"/Volumes/isomlab_1/2023.11.12.7tmps_results-filtered_hits.tar.gz"
# -------------------------------------------------------------------------

# UniProt regex (fairly strict, covers most cases):
# - 6-char format: leading [OPQ] or [A-NR-Z], then digit, then 3 alnum, then digit
_UNIPROT6 = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])'
# - 10-char format: starts with A0 (e.g., A0A…), total length 10
_UNIPROT10 = r'(?:A0[A-Z0-9]{7})'
_UNIPROT_RX = re.compile(rf'\b(?:{_UNIPROT6}|{_UNIPROT10})\b', re.IGNORECASE)

# ------------------------- Common helpers -------------------------

def _open_maybe_gz(path: str, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")

def _first_token_after_gt(header: str) -> str:
    h = header.strip()
    if not h.startswith(">"):
        return ""
    return h[1:].split(None, 1)[0] if len(h) > 1 else ""

def _uniprot_from_first_token_pipe2(tok: str) -> str:
    parts = tok.split("|")
    if len(parts) >= 2 and parts[1]:
        return parts[1].upper()
    return ""

def _uniprot_from_regex(header: str) -> str:
    m = _UNIPROT_RX.search(header)
    return m.group(0).upper() if m else ""

def _extract_uniprot(header: str) -> str:
    """Prefer 2nd '|' token from the first header token; fallback to regex anywhere in header."""
    tok = _first_token_after_gt(header)
    acc = _uniprot_from_first_token_pipe2(tok)
    if acc:
        return acc
    return _uniprot_from_regex(header)

# ------------------------- Filtering step -------------------------

def header_matches_uniprot_ids(header: str, wanted_ids_upper: Set[str]) -> bool:
    acc = _extract_uniprot(header)
    return acc in wanted_ids_upper

def filter_fasta_by_uniprot_ids(
    in_fasta: str,
    out_fasta: str,
    wanted_ids: Iterable[str],
) -> int:
    wanted_upper = {w.strip().upper() for w in wanted_ids if w and w.strip()}
    written = 0
    with _open_maybe_gz(in_fasta, "rt") as fin, open(out_fasta, "w", encoding="utf-8") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                keep = header_matches_uniprot_ids(line, wanted_upper)
                if keep:
                    fout.write(line)
                    written += 1
            else:
                if keep:
                    fout.write(line)
    return written

# ------------------------- Pairing step -------------------------

def _read_fasta_to_dict(path: str) -> Dict[str, Tuple[str, str]]:
    seqs: Dict[str, Tuple[str, str]] = {}
    with _open_maybe_gz(path, "rt") as fin:
        header = None
        chunks: List[str] = []
        for line in fin:
            if line.startswith(">"):
                if header is not None:
                    acc = _extract_uniprot(header)
                    if acc:
                        seq = "".join(chunks)
                        seqs[acc] = (header if header.endswith("\n") else header + "\n",
                                     seq if seq.endswith("\n") else seq + "\n")
                header = line.rstrip("\n")
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            acc = _extract_uniprot(header)
            if acc:
                seq = "".join(chunks)
                seqs[acc] = (header if header.endswith("\n") else header + "\n",
                             seq if seq.endswith("\n") else seq + "\n")
    return seqs

def _sanitize_for_filename(s: str) -> str:
    return "".join(ch if (ch.isalnum() or ch in "-_.") else "_" for ch in s)

def write_pairwise_fastas(
    filtered_fasta_path: str,
    uniprot_id_order: Iterable[str],
    out_dir: str,
    prefix: str = "pair",
    zero_pad: int = 6,
) -> tuple[int, list, list]:
    os.makedirs(out_dir, exist_ok=True)

    seqs = _read_fasta_to_dict(filtered_fasta_path)
    seq_keys_upper = set(seqs.keys())

    order = [uid.strip().upper() for uid in uniprot_id_order if uid and uid.strip()]
    missing = sorted({uid for uid in order if uid not in seq_keys_upper})

    written = 0
    skipped_pairs = []

    for i in range(len(order) - 1):
        a, b = order[i], order[i + 1]
        if a not in seqs or b not in seqs:
            skipped_pairs.append((a, b))
            continue

        ha, sa = seqs[a]
        hb, sb = seqs[b]

        fname = f"{prefix}_{i:0{zero_pad}d}_{_sanitize_for_filename(a)}__{_sanitize_for_filename(b)}.fa"
        out_path = os.path.join(out_dir, fname)

        with open(out_path, "w", encoding="utf-8") as fout:
            fout.write(ha); fout.write(sa)
            fout.write(hb); fout.write(sb)

        written += 1

    return written, missing, skipped_pairs

# ------------------------- ID list helpers & orchestrator -------------------------

def load_uniprot_ids_from_text(path: str) -> list[str]:
    data = Path(path).read_text(encoding="utf-8")
    items = [tok.strip() for tok in data.replace(",", "\n").splitlines()]
    return [x for x in items if x]

# ===== Alignment step (Clustal Omega) =====

"""
clustalo_pairwise_step.py
-------------------------
Drop-in step for your v2 pipeline:
  - Finds all two-sequence pair FASTA files in OUT_PAIRS_DIR
  - Runs Clustal Omega to create .aln files (Clustal format)
  - Parses alignments to compute identity/similarity metrics
  - Writes a per-pair "matches" text file mapping residue indices
  - Writes a consolidated CSV of identity metrics

Requirements:
  - Clustal Omega installed (executable name "clustalo" in PATH, or pass full path)
  - Biopython installed (AlignIO) for reading the .aln (pip install biopython)
  - No pandas dependency; CSV is written with the standard library.
"""

import csv
import os
import re
import subprocess
from pathlib import Path
from typing import Iterable, List, Tuple, Dict, Optional

# --- Clustal-like similarity symbol computation (approximate) ---
def _clustal_like_symbol(a: str, b: str) -> str:
    """Return '*', ':', '.', or ' ' for a pair of residues a,b.
       - '*' exact match (same letter, not gap, not X)
       - ':' strong similarity (both in one of strong sets)
       - '.' weak similarity (both in one of weak sets)
       - ' ' otherwise or if any gap
       Note: We do NOT exclude 'X' here; exclusion is handled in metrics.
    """
    a = a.upper()
    b = b.upper()
    if a == '-' or b == '-':
        return ' '
    if a == b and a != 'X':
        return '*'
    # Strong groups (common Clustal groupings)
    strong_groups = [
        set("STA"),
        set("NEQK"),
        set("NHQK"),
        set("NDEQ"),
        set("QHRK"),
        set("MILV"),
        set("MILF"),
        set("FYW"),
    ]
    for grp in strong_groups:
        if a in grp and b in grp:
            return ':'
    # Weak groups (approximate)
    weak_groups = [
        set("CSA"),
        set("ATV"),
        set("SAG"),
        set("STNK"),
        set("STPA"),
        set("SGND"),
        set("SNDEQK"),
        set("NDEQHK"),
        set("NEQHRK"),
        set("FVLIM"),
        set("HFY"),
    ]
    for grp in weak_groups:
        if a in grp and b in grp:
            return '.'
    return ' '

def _extract_match_lines_from_clustal(aln_path: str) -> str:
    """
    Extract the Clustal consensus symbols (* : . or space) with correct column alignment.

    We parse the file by blocks. For each block:
      - Identify the first sequence line and locate the slice that contains residues
        using a regex for [A-Za-z-]+. This defines (seq_start, seq_end) for the block.
      - Locate the consensus line (made only of spaces and '*:.'), and slice the same
        (seq_start:seq_end) region so spaces are preserved for non-matches.
    This avoids stripping trailing spaces, which would otherwise shift symbols.
    """
    blocks = []
    current = []
    with open(aln_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.strip() == "":
                if current:
                    blocks.append(current)
                    current = []
            else:
                current.append(line.rstrip("\n"))
        if current:
            blocks.append(current)

    symbols = []
    seq_re = re.compile(r"([A-Za-z\-]+)")
    for blk in blocks:
        # Find a sequence line to define slice
        seq_start = seq_end = None
        for ln in blk:
            if not ln or ln.startswith("CLUSTAL") or ln.startswith("MUSCLE") or ln.lstrip().startswith("*"):
                continue
            # Skip consensus-like lines early
            if set(ln.strip()) <= set("*:."):
                continue
            m = seq_re.search(ln)
            if m:
                seq_start, seq_end = m.span(1)
                break
        if seq_start is None:
            # No sequence slice found; skip block
            continue

        # Find consensus line (only spaces and *:. once trimmed of trailing spaces)
        cons = None
        for ln in blk:
            s = ln.rstrip()  # keep leading spaces
            if s and set(s.strip()) <= set("*:."):
                cons = ln  # original with leading spaces
        if cons is None:
            # No explicit consensus line; assume no symbols (spaces) for this block width
            width = seq_end - seq_start
            symbols.append(" " * width)
        else:
            # Append the exact symbol slice
            # If line shorter than needed, pad with spaces
            if len(cons) < seq_end:
                cons = cons + " " * (seq_end - len(cons))
            symbols.append(cons[seq_start:seq_end])

    return "".join(symbols)

# Optional Biopython import
try:
    from Bio import AlignIO
    _BIOPYTHON_OK = True
except Exception:
    _BIOPYTHON_OK = False
    AlignIO = None

def _natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\\d+)', s)]

def _list_pair_fastas(out_pairs_dir: str) -> List[str]:
    p = Path(out_pairs_dir)
    files = [str(x) for x in p.glob("*.fa")]
    files.sort(key=_natural_key)
    return files

def _ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)

def _run_clustalo(infile: str, outfile: str, clustalo_exe: str = "clustalo") -> None:
    """
    Run Clustal Omega to align a pair FASTA and write Clustal-format output.

    Parameters
    ----------
    infile : str
        Path to the input FASTA with exactly two sequences.
    outfile : str
        Path to write the Clustal (.aln) output.
    clustalo_exe : str
        Executable name or full path to Clustal Omega.
        macOS Homebrew example: '/opt/homebrew/bin/clustalo'
        Windows example (if in PATH): 'clustalo.exe'
    """
    # Ensure output directory exists
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)

    # Build command (use short options -i/-o to avoid ambiguity)
    cmd = [
        clustalo_exe,
        "-i", infile,
        "-o", outfile,
        "--seqtype=Protein",
        "--outfmt=clu",
        "--force",
    ]

    try:
        # Capture output for diagnostics; raise if non-zero exit
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        # Optional: print(res.stderr)  # Clustal writes progress to stderr
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"Clustal Omega not found: {clustalo_exe}\n"
            "On macOS with Homebrew: brew install clustal-omega\n"
            "Or pass a full path to clustalo_exe."
        ) from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"clustalo failed for {infile} with exit code {e.returncode}:\n{e.stderr}"
        ) from e


from typing import Tuple, Dict

def _parse_alignment_metrics(aln_path: str) -> Tuple[Dict[str, float], str, str, str, str, str]:
    """
    Parse a two-sequence Clustal alignment and compute identity/similarity metrics.

    IMPORTANT:
      Columns where either residue is 'X' (case-insensitive) are excluded from the
      match counters (*, :, .). Denominators (alignment length / effective length /
      seq1 length) are unchanged.
    """
    if not _BIOPYTHON_OK or AlignIO is None:
        raise RuntimeError("Biopython is required for alignment parsing. Install with: pip install biopython")

    alignment = AlignIO.read(aln_path, "clustal")
    if len(alignment) < 2:
        raise ValueError(f"Alignment has <2 sequences: {aln_path}")

    seq1_id, seq2_id = alignment[0].id, alignment[1].id
    seq1_str, seq2_str = str(alignment[0].seq), str(alignment[1].seq)

    # Use sequence length as the alignment length; Clustal may omit a consensus line in the final block.
    if len(seq1_str) != len(seq2_str):
        raise ValueError(f"Aligned sequences have different lengths: {len(seq1_str)} vs {len(seq2_str)}")

    # Denominators
    effective_length = 0
    for a, b in zip(seq1_str, seq2_str):
        if a in ["-", "X"]:
            continue
        if b in ["-", "X"]:
            continue
        effective_length += 1

    seq1_length = sum(1 for a in seq1_str if a != "-")
    seq2_length = sum(1 for a in seq2_str if a != "-")
    # Use shorter sequence for normalization
    seq_norm_length = 0
    if seq1_length >= seq2_length:
        seq_norm_length = seq2_length
    else:
        seq_norm_length = seq1_length

    # Tally matches and also build a full-length consensus line (computed, not read from file)
    full_matches = 0
    strong_similar = 0
    weak_similar = 0
    computed_match_syms = []

    for a, b in zip(seq1_str, seq2_str):
        sym = _clustal_like_symbol(a, b)  # '*' / ':' / '.' / ' '
        computed_match_syms.append(sym)

        # Skip X columns for counting (but still emit a symbol in the returned match line)
        if a.upper() == "X" or b.upper() == "X":
            continue

        if sym == "*":
            full_matches += 1
        elif sym == ":":
            strong_similar += 1
        elif sym == ".":
            weak_similar += 1

    match_lines = "".join(computed_match_syms)

    def pct(n, d): return round((n / d) * 100, 1) if d else 0.0

    metrics = {
        "Full Matches (*)": full_matches,
        "Strong Similar (:)": strong_similar,
        "Weak Similar (.)": weak_similar,

        # Ratios vs effective length (not both gaps)
        "Raw Identity (Effective) (%)": pct(full_matches, effective_length),
        "Strong Identity (Effective) (%)": pct(strong_similar, effective_length),
        "Weak Identity (Effective) (%)": pct(weak_similar, effective_length),

        # Ratios vs seq1 length
        "Raw Identity (Seq1 Norm) (%)": pct(full_matches, seq_norm_length),
        "Strong Identity (Seq1 Norm) (%)": pct(strong_similar, seq_norm_length),
        "Weak Identity (Seq1 Norm) (%)": pct(weak_similar, seq_norm_length),

        # Alignment and sequence lengths
        "Alignment Length (Effective)":effective_length,
        "Seq1 Sequence Length":seq1_length,
        "Seq2 Sequence Length":seq2_length,

    }

    metrics["Total Identity (Effective)"] = round(
        metrics["Raw Identity (Effective) (%)"]
        + metrics["Strong Identity (Effective) (%)"]
        + metrics["Weak Identity (Effective) (%)"], 1
    )
    metrics["Total Identity (Seq1 Norm)"] = round(
        metrics["Raw Identity (Seq1 Norm) (%)"]
        + metrics["Strong Identity (Seq1 Norm) (%)"]
        + metrics["Weak Identity (Seq1 Norm) (%)"], 1
    )

    return metrics, seq1_id, seq2_id, seq1_str, seq2_str, match_lines



def parse_clustalo_output_to_matches(
    seq1_str: str,
    seq2_str: str,
    match_lines: str,
    seq1_id: str,
    seq2_id: str,
) -> str:
    """
    Build a tabular mapping of aligned positions, using **PDB residue numbering**
    for each sequence when a matching PDB can be found; otherwise fall back to
    simple 1..N indexing.

    Columns:
      1) seq1_resnum_or_index ('' if gap)
      2) seq1_residue ('-' for gap)
      3) seq2_resnum_or_index ('' if gap)
      4) seq2_residue ('-' for gap)
      5) clustal_symbol (* : . or space)
    """
    # --- helpers (keep local for minimal diff) ---
    def _find_pdb_for_uid(uid: str) -> str | None:
        """Look in PDB_DIR for *uid*-trim.pdb. If not found and PDB_TARBALL exists,
        extract first matching member to OUT_PAIRS_DIR/pdb_tmp and return its path.
        """
        try:
            root = Path(PDB_DIR) if 'PDB_DIR' in globals() else None
        except Exception:
            root = None
        if root and str(root) and root.exists():
            pat = f"*{uid}*-trim.pdb"
            hits = sorted(root.rglob(pat))
            if hits:
                return str(hits[0].resolve())

        # Tarball fallback (stream & extract first match)
        try:
            tgz = Path(PDB_TARBALL) if 'PDB_TARBALL' in globals() else None
        except Exception:
            tgz = None
        if tgz and str(tgz) and tgz.exists():
            try:
                import tarfile
                out_dir = Path(OUT_PAIRS_DIR) / "pdb_tmp"
                out_dir.mkdir(parents=True, exist_ok=True)
                with tarfile.open(str(tgz), "r:gz") as tf:
                    for m in tf:
                        name = m.name.lower()
                        if not name.endswith("-trim.pdb"):
                            continue
                        if uid.lower() in name:
                            dest = out_dir / Path(m.name).name
                            # Extract only this member
                            with tf.extractfile(m) as f_in, open(dest, "wb") as f_out:
                                if f_in is not None:
                                    f_out.write(f_in.read())
                            return str(dest.resolve())
            except Exception as e:
                print(f"[WARN] Tarball scan failed for {uid}: {e}")

        return None

    def _parse_pdb_residue_series(pdb_path: str) -> list[str]:
        """Return ordered residue identifiers like ['1','2','42A',...] from ATOM lines.
        Uses first chain encountered.
        """
        res = []
        seen = set()
        chain_chosen = None
        try:
            with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line.startswith("ATOM"):
                        continue
                    chain = line[21].strip()
                    if chain_chosen is None:
                        chain_chosen = chain
                    if chain_chosen != chain:
                        continue
                    resseq = line[22:26].strip()
                    icode  = line[26].strip()
                    if not resseq:
                        continue
                    key = (resseq, icode)
                    if key in seen:
                        continue
                    seen.add(key)
                    res.append(f"{resseq}{icode}" if icode else f"{resseq}")
        except Exception as e:
            print(f"[WARN] Failed to parse PDB residues from {pdb_path}: {e}")
        return res

    def _map_alignment_to_numbers(aligned_seq: str, series: list[str] | None) -> list[str | None]:
        """Consume next residue id for each non-gap; gaps -> None. If series is None, use 1..N indexing."""
        out: list[str | None] = []
        if series is None:
            # simple 1..N indexing by non-gap letters
            idx = 0
            for ch in aligned_seq:
                if ch == '-':
                    out.append(None)
                else:
                    idx += 1
                    out.append(str(idx))
            return out
        i = 0
        n = len(series)
        for ch in aligned_seq:
            if ch == '-':
                out.append(None)
            else:
                if i < n:
                    out.append(series[i])
                    i += 1
                else:
                    out.append(None)  # ran out of PDB residues
        return out

    # Try to find PDBs for each UniProt id
    pdb1 = _find_pdb_for_uid(seq1_id)
    pdb2 = _find_pdb_for_uid(seq2_id)
    series1 = _parse_pdb_residue_series(pdb1) if pdb1 else None
    series2 = _parse_pdb_residue_series(pdb2) if pdb2 else None

    map1 = _map_alignment_to_numbers(seq1_str, series1)
    map2 = _map_alignment_to_numbers(seq2_str, series2)

    # Build rows
    out_lines = []
    for a, b, r1, r2, sym in zip(seq1_str, seq2_str, map1, map2, match_lines):
        col1 = "" if a == "-" else (r1 or "")
        col2 = a
        col3 = "" if b == "-" else (r2 or "")
        col4 = b
        out_lines.append(f"{col1}\t{col2}\t{col3}\t{col4}\t{sym}")

    header = [
        f"# Column 1 and 2 map to {seq1_id} ({'PDB residue numbers' if series1 else '1..N indices'})",
        f"# Column 3 and 4 map to {seq2_id} ({'PDB residue numbers' if series2 else '1..N indices'})",
        f"# Column 5 indicates Clustal match: exact(*), strong(:), weak(.), space = none",
    ]
    return "\n".join(header) + "\n" + "\n".join(out_lines) + "\n"

def run_clustalo_on_pairs(
    out_pairs_dir: str,
    clustalo_exe: str = "clustalo",
    matches_dir: Optional[str] = None,
    csv_basepath: Optional[str] = None,
) -> str:
    """
    Run on all *.fa in out_pairs_dir (sorted naturally).
    - Writes *_aligned.aln next to each .fa
    - Writes a matches text file per pair into matches_dir (default: out_pairs_dir / 'matches')
    - Writes a consolidated CSV (csv_basepath + '_traceback_identity.csv').
      If csv_basepath is None, uses out_pairs_dir.
    Returns the CSV path.
    """
    pair_files = _list_pair_fastas(out_pairs_dir)
    if not pair_files:
        raise FileNotFoundError(f"No .fa pairs found in {out_pairs_dir}")

    matches_dir = matches_dir or str(Path(out_pairs_dir) / "matches")
    _ensure_dir(matches_dir)

    csv_basepath = csv_basepath or out_pairs_dir
    csv_path = str(Path(csv_basepath).with_suffix("")) + "_traceback_identity.csv"

    identity_rows = []
    for i, fasta_path in enumerate(reversed(pair_files), start=0):
        aln_path = fasta_path.replace(".fa", "_aligned.aln")
        # Align
        print(f"[clustalo] {i}/{len(pair_files)}: {os.path.basename(fasta_path)}")
        _run_clustalo(fasta_path, aln_path, clustalo_exe=clustalo_exe)

        # Parse metrics
        metrics, s1_id, s2_id, s1_str, s2_str, match_lines = _parse_alignment_metrics(aln_path)
        s1_id = s1_id.split("|")[1]
        s2_id = s2_id.split("|")[1]

        # Build matches report
        matches_filename = f"{i}_{s1_id}*{s2_id}_matches.txt"
        matches_path = str(Path(matches_dir) / matches_filename)
        with open(matches_path, "w", encoding="utf-8") as fh:
            fh.write(parse_clustalo_output_to_matches(s1_str, s2_str, match_lines, s1_id, s2_id))

        # Consolidate metrics row
        row = {
            "Traceback Depth": i,
            "Sequence 1": s1_id,
            "Sequence 2": s2_id,
        }
        row.update(metrics)
        identity_rows.append(row) 

    # Write CSV
    field_order = [
        "Traceback Depth", "Sequence 1", "Sequence 2",
        "Full Matches (*)", "Strong Similar (:)", "Weak Similar (.)",
        "Total Identity (Effective)",
        "Raw Identity (Effective) (%)",
        "Strong Identity (Effective) (%)",
        "Weak Identity (Effective) (%)",
        "Total Identity (Seq1 Norm)",
        "Raw Identity (Seq1 Norm) (%)", "Strong Identity (Seq1 Norm) (%)", "Weak Identity (Seq1 Norm) (%)",
        "Alignment Length (Effective)", "Seq1 Sequence Length", "Seq2 Sequence Length",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=field_order)
        w.writeheader()
        for row in identity_rows:
            w.writerow(row)

    print(f"[DONE] Consolidated identity report saved to {csv_path}")
    return csv_path



def build_uniprot_to_matches_map(matches_dir: str) -> Dict[str, str]:
    """
    Build a dict mapping **both** UniProt IDs present in each *_matches.txt file name
    to that file's absolute path. This ensures the final pair contributes two mappings
    (one per PDB) when both have corresponding structures.
    """
    out: Dict[str, str] = {}
    p = Path(matches_dir)
    if not p.exists():
        print(f"[WARN] matches_dir does not exist: {matches_dir}")
        return out

    # Two UniProt patterns (6-char and 10-char, uppercase)
    _UNIPROT6 = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])'
    _UNIPROT10 = r'(?:[A-NR-Z0-9]{4}[A-Z0-9]{6})'  # coarse but works for A0A...

    uid_re = re.compile(rf'({_UNIPROT10}|{_UNIPROT6})')

    for f in sorted(p.glob("*_matches.txt")):
        name = f.name
        # Extract in-order appearances to preserve pair ordering (seq1, seq2)
        uids = uid_re.findall(name.upper())
        uids_unique = []
        for u in uids:
            if u not in uids_unique:
                uids_unique.append(u)
        if not uids_unique:
            # Last resort: try to read first line to recover IDs (optional)
            pass

        # Map both IDs (if present) to the same match file
        abs_path = str(f.resolve())
        for u in uids_unique[:2]:
            out[u] = abs_path

    return out

def open_trimmed_pdbs_from_tar(
    tar_gz_path: str,
    allowed_uniprot_ids: Iterable[str],
    out_tmp_dir: str
) -> Dict[str, Dict]:
    """
    Stream a .tar.gz of PDB files, selecting only files whose names end with '-trim.pdb'
    and that contain a UniProt ID in `allowed_uniprot_ids`. For performance with very
    large archives (e.g., 70GB), this function:
      - Uses streaming tar modes (no in-memory member lists).
      - Optionally uses pigz for multi-threaded gzip decompression if available.
      - Copies matched members to a small temp dir and opens them with PDBfile.
    """
    out: Dict[str, Dict] = {}
    allowed = {uid.upper() for uid in allowed_uniprot_ids if uid}
    out_dir = Path(out_tmp_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not tar_gz_path:
        print("[WARN] No tar.gz path specified; skipping PDB open step.")
        return out
    tgz = Path(tar_gz_path)
    if not tgz.exists():
        print(f"[WARN] PDB tarball not found: {tar_gz_path}")
        return out
    if not _PDBFILE_AVAILABLE:
        print("[WARN] pdbFile.PDBfile is not importable; skipping PDB open step.")
        return out

    import tarfile as _tarfile
    import io, subprocess, shutil

    # Prefer pigz for faster gzip decompression (multi-threaded), if available.
    pigz = shutil.which("pigz")
    if pigz:
        # pigz -dc file.tgz | tarfile.open(fileobj=proc.stdout, mode="r|")
        proc = subprocess.Popen([pigz, "-dc", str(tgz)], stdout=subprocess.PIPE)
        fileobj = proc.stdout
        tar_mode = "r|"
    else:
        # Fallback: Python's single-threaded gzip; still streaming
        proc = None
        fileobj = None
        tar_mode = "r|gz"

    try:
        with _tarfile.open(tar_gz_path if fileobj is None else None, mode=tar_mode, fileobj=fileobj) as tf:
            for m in tf:  # streaming iteration (no getmembers)
                if not m.isfile():
                    continue
                name = m.name
                # Require '-trim.pdb' suffix (case-insensitive)
                if not name.lower().endswith("-trim.pdb"):
                    continue

                # Fast allowlist check: substring presence first
                uid_match = None
                for uid in allowed:
                    if uid in name:
                        uid_match = uid
                        break
                if not uid_match:
                    # Fallback: regex scan + allowlist check
                    m_u = _UNIPROT_RX.search(name)
                    if m_u:
                        cand = m_u.group(0).upper()
                        if cand in allowed:
                            uid_match = cand
                if not uid_match:
                    continue

                # Extract just this member as a stream and write in chunks
                fobj = tf.extractfile(m)
                if fobj is None:
                    continue
                fname = Path(name).name
                tmp_path = out_dir / fname
                with open(tmp_path, "wb") as fout:
                    shutil.copyfileobj(fobj, fout, length=1024 * 1024)  # 1MB chunks

                # Open with PDBFile
                try:
                    pdb = _PDBfile(pdbFilePath=str(tmp_path.parent) + str(Path('/')),
                                   pdbFileName=tmp_path.name,
                                   pdbFileAsString="",
                                   twoCharacterChain=0,
                                   zip_status=0)
                    n_atoms = len(getattr(pdb, "atoms", {})) + len(getattr(pdb, "hetatoms", {}))
                    n_res = len(getattr(pdb, "residues", {})) + len(getattr(pdb, "het_residues", {}))
                    chains = list(getattr(pdb, "chains", []))
                    out[uid_match] = {
                        "pdb_path": str(tmp_path.resolve()),
                        "pdb_code": getattr(pdb, "pdbCode", Path(tmp_path).stem),
                        "n_atoms": n_atoms,
                        "n_residues": n_res,
                        "chains": chains,
                        "tar_member": name,
                    }
                except Exception as e:
                    print(f"[WARN] Failed to open PDB for {uid_match} ({name}): {e}")
                    out[uid_match] = {
                        "pdb_path": str(tmp_path.resolve()),
                        "pdb_code": Path(tmp_path).stem,
                        "error": str(e),
                        "tar_member": name,
                    }
    finally:
        if pigz and proc:
            try:
                proc.stdout.close()
            except Exception:
                pass
            proc.wait()

    print(out)
    return out


def open_trimmed_pdbs_from_dir(
    pdb_dir: str,
    allowed_uniprot_ids: Iterable[str],
) -> Dict[str, Dict]:
    """
    Scan a directory tree of PDB files and open only files whose names end with '-trim.pdb'
    and whose names contain a UniProt ID in `allowed_uniprot_ids`. Uses pdbFile.PDBfile.
    Returns a dict keyed by UniProt ID with a small summary.
    """
    out: Dict[str, Dict] = {}
    if not pdb_dir:
        print("[WARN] No PDB_DIR specified; skipping PDB dir step.")
        return out
    root = Path(pdb_dir)
    if not root.exists():
        print(f"[WARN] PDB_DIR not found: {pdb_dir}")
        return out
    if not _PDBFILE_AVAILABLE:
        print("[WARN] pdbFile.PDBfile is not importable; skipping PDB open step.")
        return out

    allowed = {uid.upper() for uid in allowed_uniprot_ids if uid}
    for dirpath, _dnames, fnames in os.walk(root):
        for fname in fnames:
            name_lower = fname.lower()
            if not name_lower.endswith("-trim.pdb"):
                continue
            # Prefer direct substring match against allowed IDs
            uid_match = None
            for uid in allowed:
                if uid in fname:
                    uid_match = uid
                    break
            if not uid_match:
                # Fallback: regex on filename
                m_u = _UNIPROT_RX.search(fname)
                if m_u:
                    cand = m_u.group(0).upper()
                    if cand in allowed:
                        uid_match = cand
            if not uid_match:
                continue

            fpath = Path(dirpath) / fname
            try:
                pdb = _PDBfile(pdbFilePath=str(fpath.parent) + os.sep,
                               pdbFileName=fpath.name,
                               pdbFileAsString="",
                               twoCharacterChain=0,
                               zip_status=0)
                n_atoms = len(getattr(pdb, "atoms", {})) + len(getattr(pdb, "hetatoms", {}))
                n_res = len(getattr(pdb, "residues", {})) + len(getattr(pdb, "het_residues", {}))
                chains = list(getattr(pdb, "chains", []))
                out[uid_match] = {
                    "pdb_path": str(fpath.resolve()),
                    "pdb_code": getattr(pdb, "pdbCode", fpath.stem),
                    "n_atoms": n_atoms,
                    "n_residues": n_res,
                    "chains": chains,
                    "tar_member": None,
                }
            except Exception as e:
                print(f"[WARN] Failed to open PDB for {uid_match} ({fpath.name}): {e}")
                out[uid_match] = {
                    "pdb_path": str(fpath.resolve()),
                    "pdb_code": fpath.stem,
                    "error": str(e),
                    "tar_member": None,
                }
    # print(out)  # print dictionary for inspection
    return out


from pathlib import Path
from typing import List, Tuple, Optional

def extract_sequence_conservation(file_path: str, file_name: str, return_position: int = 2) -> List[Tuple[int, Optional[int], str]]:
    """
    Parse a pairwise match file to extract sequence conservation information.

    Accepts either:
      - 5-column lines: pos1, aa1, pos2, aa2, symbol
      - 3-column lines: pos1, pos2, symbol

    Args:
        file_path (str): Directory containing the file (can be "" for current dir).
        file_name (str): Match filename.
        return_position (int): Which position to treat as primary (1 or 2). Default=2.

    Returns:
        List[Tuple[int, Optional[int], str]]:
            A list of (selected_position, other_position, symbol), where
            selected_position corresponds to return_position (1 or 2), and
            other_position is the remaining one (may be None if missing).
    """
    # Normalize which position to return
    rp = 1 if return_position == 1 else 2

    # Build path robustly (works even if file_path is empty or lacks a trailing slash)
    pairwise_file = str((Path(file_path) / file_name) if file_path else Path(file_name))

    results: List[Tuple[int, Optional[int], str]] = []
    try:
        with open(pairwise_file, "r", encoding="utf-8") as match_file:
            for line in match_file:
                # Skip comments or blank lines
                if not line.strip() or line.lstrip().startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t")
                pos1: Optional[int] = None
                pos2: Optional[int] = None
                symbol: Optional[str] = None

                if len(parts) >= 5:
                    # 5-col format: pos1, aa1, pos2, aa2, symbol
                    try:
                        pos1 = int(parts[0]) if parts[0] not in ("-", "") else None
                        pos2 = int(parts[2]) if parts[2] not in ("-", "") else None
                        symbol = parts[4]
                    except ValueError:
                        continue
                elif len(parts) >= 3:
                    # 3-col format: pos1, pos2, symbol
                    try:
                        pos1 = int(parts[0]) if parts[0] not in ("-", "") else None
                        pos2 = int(parts[1]) if parts[1] not in ("-", "") else None
                        symbol = parts[2]
                    except ValueError:
                        continue
                else:
                    continue

                if symbol is None:
                    continue

                # Choose which position is "selected" and which is "other"
                selected = pos1 if rp == 1 else pos2
                other    = pos2 if rp == 1 else pos1

                # Require the selected position to exist; skip rows where it's missing
                if selected is None:
                    continue

                results.append((selected, other, symbol))
    except FileNotFoundError:
        print(f"[WARN] Match file not found: {pairwise_file}")

    return results



def map_sequence_conservation(pdb_file_path: str, pdb_file_name: str, out_path : str, position_matches):
    """
    Opens the given PDB using pdbFile.PDBfile and maps conservation symbols onto
    atom temperature factors for residues in the first sequence.
    Returns the PDBfile object with updated temperature factors.
    """
    if not _PDBFILE_AVAILABLE:
        print("[WARN] pdbFile.PDBfile not available; skipping map.")
        return None

    # Import here to avoid top-level import errors if unavailable
    import os as _os
    try:
        pdb = _PDBfile(pdbFilePath=str(pdb_file_path) + (_os.sep if not str(pdb_file_path).endswith((_os.sep, '/')) else ""),
                       pdbFileName=pdb_file_name,
                       pdbFileAsString="",
                       twoCharacterChain=0,
                       zip_status=0)
    except Exception as e:
        print(f"[WARN] Failed to open PDB: {pdb_file_path}{_os.sep}{pdb_file_name}: {e}")
        return None

    # Build a quick index from residue number -> residue object (first CA found wins)
    # PDBfile exposes residues as a dict; each residue has .num and .atoms
    num_to_res = {}
    for key in pdb.residues:
        try:
            res = pdb.residues[key]
            # print(key)
            # print(res)
            num_to_res.setdefault(res.num, res)
        except Exception:
            continue

    # for x in position_matches:
    #     print(x)


    # print(pdb_file_name)
    # raise SystemExit

    # Set values for matching residues
    matched_nums = set()
    for pos1, _pos2, symbol in position_matches:
        res = num_to_res.get(pos1)
        # print(pos1, symbol, res)
        if not res:
            continue
        # print(pos1, symbol)
        matched_nums.add(pos1)
        for atom_key, atom in getattr(res, "atoms", {}).items():
            try:
                if symbol == "*":
                    atom.temperature_factor = 100
                elif symbol == ":":
                    atom.temperature_factor = 50
                elif symbol == ".":
                    atom.temperature_factor = 25
                else:
                    # non-matching or undefined symbol
                    atom.temperature_factor = 0
            except Exception:
                pass

    # print(pdb_file_name)
    # raise SystemExit

    # Define non-matching residues as 0
    for key in pdb.residues:
        try:
            res = pdb.residues[key]
            if res.num not in matched_nums:
                for atom_key, atom in getattr(res, "atoms", {}).items():
                    atom.temperature_factor = 0
        except Exception:
            continue

    # Ensure output directory ends with path separator and exists
    if out_path:
        if not out_path.endswith((_os.sep, '/')):
            out_path = out_path + _os.sep
        try:
            Path(out_path).mkdir(parents=True, exist_ok=True)
        except Exception:
            pass

    # Write the revised PDB file
    out_file = out_path + pdb_file_name.split(".pdb")[0] + "_mc.pdb"
    with open(out_file, "w", encoding="utf-8") as f:
        for key in pdb.residues:
            f.write(str(pdb.residues[key]))

    return pdb

def run_pipeline(
    in_fasta: str,
    out_filtered_fasta: str,
    wanted_ids_path: str,
    out_pairs_dir: str,
    pair_filename_prefix: str = "pair",
    zero_pad: int = 6,
) -> None:
    print("[STEP 1/2] Loading UniProt ID order...")
    wanted_ids = load_uniprot_ids_from_text(wanted_ids_path)
    print(f"[INFO] IDs in list: {len(wanted_ids):,}")

    print("[STEP 1/2] Filtering FASTA (pipe-2nd / regex fallback)...")
    n = filter_fasta_by_uniprot_ids(in_fasta, out_filtered_fasta, wanted_ids)
    print(f"[INFO] Wrote {n:,} matching records to: {out_filtered_fasta}")

    print("[STEP 2/2] Writing pairwise FASTA files...")
    num_pairs, missing_ids, skipped = write_pairwise_fastas(
        out_filtered_fasta, wanted_ids, out_pairs_dir,
        prefix=pair_filename_prefix, zero_pad=zero_pad
    )
    print(f"[INFO] Wrote {num_pairs:,} pair files → {out_pairs_dir}")

    if missing_ids:
        print(f"[WARN] {len(missing_ids):,} IDs from the list were not present in the filtered FASTA.")
        if len(missing_ids) <= 20:
            print("       Missing:", ", ".join(missing_ids))
        else:
            print("       (showing first 20)", ", ".join(missing_ids[:20]))

    if skipped:
        print(f"[WARN] Skipped {len(skipped):,} pairs due to missing sequences.")
        if len(skipped) <= 10:
            print("       Skipped pairs:", skipped)
        else:
            print("       (showing first 10)", skipped[:10])

    print("[DONE] Pipeline complete.")

    # Alignment step: run Clustal Omega on generated pairs
    csv_path = run_clustalo_on_pairs(
        out_pairs_dir=OUT_PAIRS_DIR,
        clustalo_exe="clustalo",
        matches_dir=str(Path(OUT_PAIRS_DIR) / "matches"),
        csv_basepath=str(Path(OUT_PAIRS_DIR) / "matches")
    )
    print(f"[DONE] Alignment metrics CSV: {csv_path}")
    matches_dir_path = Path(OUT_PAIRS_DIR) / "matches"
    _uid_to_match = build_uniprot_to_matches_map(str(matches_dir_path))
    try:
        csv_path = run_clustalo_on_pairs(
            out_pairs_dir=OUT_PAIRS_DIR,
            clustalo_exe="clustalo",
            matches_dir=str(Path(OUT_PAIRS_DIR) / "matches"),
            csv_basepath=str(Path(OUT_PAIRS_DIR) / "matches")
        )
        print(f"[DONE] Alignment metrics CSV: {csv_path}")
        matches_dir_path = Path(OUT_PAIRS_DIR) / "matches"
        _uid_to_match = build_uniprot_to_matches_map(str(matches_dir_path))
        # === PDB step ===
        try:
            matches_dir_path = Path(OUT_PAIRS_DIR) / "matches"
            _uid_to_match = build_uniprot_to_matches_map(str(matches_dir_path))
            print(f"[INFO] UniProt→matches entries: {len(_uid_to_match)}")
            # Prefer directory scan if PDB_DIR is set and exists; else fall back to tarball
            if PDB_DIR and Path(PDB_DIR).exists():
                _uid_to_pdbinfo = open_trimmed_pdbs_from_dir(PDB_DIR, _uid_to_match.keys())
            elif PDB_TARBALL and Path(PDB_TARBALL).exists():
                pdb_tmp_dir = str(Path(OUT_PAIRS_DIR) / "pdb_tmp")
                _uid_to_pdbinfo = open_trimmed_pdbs_from_tar(PDB_TARBALL, _uid_to_match.keys(), pdb_tmp_dir)
            else:
                print("[WARN] Neither PDB_DIR nor PDB_TARBALL is available; skipping PDB step.")
                _uid_to_pdbinfo = {}
            print(f"[INFO] Opened PDBs: {len(_uid_to_pdbinfo)}")

        except Exception as e:
            print("[ERROR] PDB step failed:", e)
        # === Map sequence conservation onto PDBs ===
        try:
            common_ids = sorted(set(_uid_to_match.keys()) & set(_uid_to_pdbinfo.keys()))
            print(f"[INFO] Common UniProt IDs (matches ∩ PDBs): {len(common_ids)}")
            annot_dir = Path(OUT_PAIRS_DIR) / "annotated_pdbs"
            annot_dir.mkdir(parents=True, exist_ok=True)
            _uid_to_pdb_annot = {}
            
            
            # Special-case: dual-map ONLY the last *_matches.txt file in order
            try:
                _all_matches = sorted(Path(matches_dir_path).glob("*_matches.txt"))
                _last_match = _all_matches[-1] if _all_matches else None
                if _last_match is not None:
                    _name = _last_match.name.upper()
                    _uid_order = re.findall(r'(?:A0A[A-Z0-9]{7}|[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])', _name)
                    # Deduplicate while preserving order
                    _uid_order_dedup = []
                    for _u in _uid_order:
                        if _u not in _uid_order_dedup:
                            _uid_order_dedup.append(_u)
                    if len(_uid_order_dedup) >= 2:
                        _first_uid, _second_uid = _uid_order_dedup[0], _uid_order_dedup[1]
                        # Only proceed for UIDs that actually have PDBs opened
                        for _u, _rp in [(_first_uid, 1), (_second_uid, 2)]:
                            if _u in _uid_to_pdbinfo and _u in _uid_to_match:
                                _mpath = Path(_uid_to_match[_u])
                                _pos_matches = extract_sequence_conservation(str(_mpath.parent) + str(Path('/')), _mpath.name, return_position=_rp)
                                _ppath = Path(_uid_to_pdbinfo[_u]["pdb_path"])
                                _pdb_obj = map_sequence_conservation(str(_ppath.parent), _ppath.name, str(annot_dir), _pos_matches)
                                if _pdb_obj is not None:
                                    _uid_to_pdb_annot[_u] = _pdb_obj
                        # Exclude these two from the regular loop to avoid duplicates
                        common_ids = {c for c in common_ids if c not in {_first_uid, _second_uid}}
            except Exception as _e:
                print("[WARN] Dual-map (last file) step skipped due to error:", _e)
            for uid in common_ids:
                # 1) Parse match file -> position matches
                match_abs = _uid_to_match[uid]
                mpath = Path(match_abs)
                
                # Determine which position in the pair this UID occupies
                # by parsing the match filename for two UniProt IDs in order.
                _uid_order = []
                try:
                    _uid_order = re.findall(r'(?:A0A[A-Z0-9]{7}|[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])', mpath.name.upper())
                except Exception:
                    _uid_order = []
                _rp = 2
                if len(_uid_order) >= 2:
                    if uid.upper() == _uid_order[0]:
                        _rp = 1
                    elif uid.upper() == _uid_order[1]:
                        _rp = 2
                # Parse matches with the correct orientation for this UID
                pos_matches = extract_sequence_conservation(str(mpath.parent) + str(Path('/')), mpath.name, return_position=_rp)


                # 2) Open PDB for uid and map conservation
                pdb_abs = _uid_to_pdbinfo[uid]["pdb_path"]
                ppath = Path(pdb_abs)
                pdb_obj = map_sequence_conservation(str(ppath.parent), ppath.name, str(annot_dir), pos_matches)
                if pdb_obj is not None:
                    _uid_to_pdb_annot[uid] = pdb_obj

            print(f"[INFO] Annotated PDBs: {len(_uid_to_pdb_annot)}")
        except Exception as e:
            print("[ERROR] Conservation mapping step failed:", e)
    except Exception as e:
        print("[ERROR] Alignment step failed:", e)


if __name__ == "__main__":
    if IN_FASTA and OUT_FILTERED_FASTA and WANTED_IDS_PATH and OUT_PAIRS_DIR:
        run_pipeline(IN_FASTA, OUT_FILTERED_FASTA, WANTED_IDS_PATH, OUT_PAIRS_DIR)
    else:
        print("Edit the CONFIG section at the top of this file, then run again.")
        print("Or import this module and call run_pipeline(in_fasta, out_filtered_fasta, wanted_ids_path, out_pairs_dir).")

