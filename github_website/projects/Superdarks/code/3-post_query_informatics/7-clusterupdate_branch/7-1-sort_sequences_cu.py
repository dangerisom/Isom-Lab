#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
split_fasta_by_domain.py
------------------------
Splits one FASTA into:
  <base>.archaea.fa
  <base>.bacteria.fa
  <base>.eukaryota.fa
  <base>.other.fa
and writes a report file:
  <base>.split_report.txt

Classification looks at the first field of the FASTA header before the first '|'
(e.g., '>Bacteria|…', '>Eukaryota|…', '>Archaea|…'), case-insensitive. Falls back
to scanning the whole header if needed. Works with .fa/.fasta and .gz variants.
"""

from pathlib import Path
import gzip

# ===== EDIT THIS to your file name (supports .fa/.fasta and .gz variants) =====
INPUT_PATH = Path("/Users/danisom/Desktop/local_informatics/2023.11.12.trim_blastp_sequences_final.fa")
# ==============================================================================

def base_stem(p: Path) -> str:
    name = p.name
    for suf in (".fa.gz", ".fasta.gz", ".fna.gz", ".faa.gz"):
        if name.lower().endswith(suf):
            return name[:-len(suf)]
    for suf in (".fa", ".fasta", ".fna", ".faa"):
        if name.lower().endswith(suf):
            return name[:-len(suf)]
    return p.stem

def open_text_maybe_gz(path: Path, mode="rt"):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")

def classify_header(header_line: str) -> str:
    """Return one of: 'archaea','bacteria','eukaryota','other'."""
    h = header_line.lstrip(">").strip()
    # Prefer first '|' field if present, else first whitespace token
    first = h.split("|", 1)[0].split(None, 1)[0].lower()

    # Normalize common variants
    if first in ("archaea", "archaeal", "archaeota", "archaeon"):
        return "archaea"
    if first in ("bacteria", "bacterial", "bacterium"):
        return "bacteria"
    if first in ("eukaryota", "eukaryotic", "eukaryote", "eukarya"):
        return "eukaryota"

    # Fallback: search anywhere in header
    low = h.lower()
    if "archaea" in low or "archaeota" in low:
        return "archaea"
    if "bacteria" in low:
        return "bacteria"
    if any(x in low for x in ("eukaryota","eukaryote","eukaryotic","eukarya")):
        return "eukaryota"
    return "other"

def main():
    in_path = Path(INPUT_PATH)
    if not in_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {in_path}")

    base = base_stem(in_path)
    out_dir = in_path.parent

    out_paths = {
        "archaea":   out_dir / f"{base}.archaea.fa",
        "bacteria":  out_dir / f"{base}.bacteria.fa",
        "eukaryota": out_dir / f"{base}.eukaryota.fa",
        "other":     out_dir / f"{base}.other.fa",
    }
    report_path = out_dir / f"{base}.split_report.txt"

    # Open all outputs up front (keep empty files so report paths are valid)
    outs = {k: open(v, "w", encoding="utf-8", newline="\n") for k, v in out_paths.items()}
    counts = {k: 0 for k in out_paths}

    try:
        with open_text_maybe_gz(in_path, "rt") as fin:
            current_key = None
            current_lines = []
            for line in fin:
                if line.startswith(">"):
                    # flush previous record
                    if current_key is not None and current_lines:
                        outs[current_key].writelines(current_lines)
                    # start new record
                    current_key = classify_header(line)
                    counts[current_key] += 1
                    current_lines = [line]
                else:
                    if current_key is None:
                        # skip until first header
                        continue
                    current_lines.append(line)
            # flush last record
            if current_key is not None and current_lines:
                outs[current_key].writelines(current_lines)
    finally:
        for f in outs.values():
            f.close()

    # Write the report with the exact look you requested
    # Two leading spaces; left-justify key to width 11, then arrow, file name, and counts.
    with open(report_path, "w", encoding="utf-8", newline="\n") as rep:
        for key in ("archaea", "bacteria", "eukaryota", "other"):
            rep.write(f"  {key:<11} -> {out_paths[key].name} ({counts[key]} seqs)\n")

    # Also print a quick summary to console
    print("Done. Wrote:")
    for key in ("archaea","bacteria","eukaryota","other"):
        print(f"  {key:<11} -> {out_paths[key].name} ({counts[key]} seqs)")
    print(f"Report: {report_path.name}")

if __name__ == "__main__":
    main()
