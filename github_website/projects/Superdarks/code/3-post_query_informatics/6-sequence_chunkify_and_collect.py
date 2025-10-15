import csv
import re
from pathlib import Path
from typing import Dict, Tuple

# Columns to pull from the TSV (keys = output names, values = TSV header names)
_REQUIRED_COLS = {
    "rank": "rank",
    "uniprot_id": "primary_accession",
    "gene": "gene",
    "protein": "protein",
    "organism": "organism",
    "superkingdom": "superkingdom",
}

# Treat these as "empty" → UC (for gene/protein only)
_EMPTY_TOKENS = {"", "na", "n/a", "none", ".", "null"}

# Match ALL unicode whitespace (space, tabs, NBSP, etc.)
_WS_RE = re.compile(r"\s+", flags=re.UNICODE)

def _normalize_name(value: str) -> str:
    """
    For gene/protein:
      - if empty/NA-like -> 'UC'
      - normalize odd spaces -> underscores
    """
    if value is None:
        return "UC"
    v = str(value).strip()
    if v.lower() in _EMPTY_TOKENS:
        return "UC"
    v = v.replace("\u00A0", " ").replace("\u2007", " ").replace("\u202F", " ")
    v = _WS_RE.sub("_", v)
    v = re.sub(r"_+", "_", v).strip("_")
    return v if v else "UC"

def _normalize_organism(value: str) -> str:
    """
    For organism:
      - keep empty as empty (no 'UC' unless you want it)
      - normalize odd spaces -> underscores
    """
    if value is None:
        return ""
    v = str(value).strip()
    if not v:
        return ""
    v = v.replace("\u00A0", " ").replace("\u2007", " ").replace("\u202F", " ")
    v = _WS_RE.sub("_", v)
    v = re.sub(r"_+", "_", v).strip("_")
    return v

def _extract_entry(row: Dict[str, str]) -> Tuple[str, Dict[str, str]]:
    """
    Convert one DictReader row to (uniprot_id, entry_dict).
    Ensures gene/protein sanitized (UC + underscores) and organism underscores.
    """
    norm = {k: (v.strip() if isinstance(v, str) else ("" if v is None else str(v)))
            for k, v in row.items()}

    uid = norm.get(_REQUIRED_COLS["uniprot_id"], "")
    if not uid:
        return "", {}  # skip rows with no UniProt ID

    entry = {
        "rank": norm.get(_REQUIRED_COLS["rank"], ""),
        "uniprot_id": uid,
        "gene": _normalize_name(norm.get(_REQUIRED_COLS["gene"], "")),
        "protein": _normalize_name(norm.get(_REQUIRED_COLS["protein"], "")),
        "organism": _normalize_organism(norm.get(_REQUIRED_COLS["organism"], "")),
        "superkingdom": norm.get(_REQUIRED_COLS["superkingdom"], ""),
    }
    return uid, entry

def build_uniprot_index(tsv_path: Path) -> Dict[str, Dict[str, str]]:
    """
    Build a dict:
      { <UniProt ID>: {
          "rank", "uniprot_id", "gene", "protein", "organism", "superkingdom"
        }, ... }
    Keeps the FIRST occurrence for duplicate UniProt IDs.
    """
    index: Dict[str, Dict[str, str]] = {}
    with tsv_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [v for v in _REQUIRED_COLS.values() if v not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"TSV is missing required columns: {missing}")

        for row in reader:
            uid, entry = _extract_entry(row)
            if not uid:
                continue
            if uid in index:
                continue  # change to: index[uid] = entry  if you prefer "last wins"
            index[uid] = entry
    return index

def extract_protein_sequence_with_explicit_gaps(pdb_file_path):
    """
    Extracts the protein sequence from a PDB file using CA atoms, preserving gaps
    by inserting '-' for every missing residue in sequential numbering.
    
    Args:
        pdb_file_path (str): Path to the PDB file.
    
    Returns:
        str: The extracted protein sequence in single-letter amino acid code, with gaps.
    """
    # Map of three-letter codes to one-letter codes
    three_to_one = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
    }

    sequence = []
    previous_residue_number = None

    with open(pdb_file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM") and " CA " in line:
                # Extract residue name (columns 17-20 in PDB format)
                residue_name = line[17:20].strip()
                # Extract residue sequence number (columns 23-26 in PDB format)
                residue_number = int(line[22:26].strip())

                # Add gaps for missing residues
                if previous_residue_number is not None and residue_number > previous_residue_number + 1:
                    gap_size = residue_number - previous_residue_number - 1
                    sequence.extend(["X"] * gap_size)

                # Map to single-letter code if valid
                if residue_name in three_to_one:
                    sequence.append(three_to_one[residue_name])

                # Update the previous residue number
                previous_residue_number = residue_number

    # Join the sequence list into a single
    return "".join(sequence)

if __name__ == "__main__":
    import os
    import math
    import subprocess
    from typing import List
    from pathlib import Path

    # ------------------- CONFIG (edit to taste) -------------------
    TSV_FILE = Path("/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/2023.11.12.7tmps-analysis_final/output_with_interpro.tsv")
    PDB_FILE_PATH = Path("/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/2023.11.12.7tmps-analysis_final/2023.11.12.7tmps_results-filtered_hits/")
    OUT_DIR = Path("/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/2023.11.12.7tmps-analysis_final")

    N_CHUNKS = 256                # total chunks to split the file list into
    JOB_NAME = "trimseq"          # LSF job family name
    QUEUE = "isomlab"             # LSF queue name

    # ----------------- helpers: chunking + I/O per chunk -----------
    def _chunkify(items, n_chunks):
        """Split list into n nearly-equal chunks; last chunk may be shorter."""
        items = list(items)
        if not items or n_chunks < 1:
            return []
        size = math.ceil(len(items) / n_chunks)
        return [items[i*size:(i+1)*size] for i in range(n_chunks) if items[i*size:(i+1)*size]]

    def _write_chunk_lists(chunks_dir: Path, chunks: List[List[str]]) -> List[Path]:
        chunks_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for i, chunk in enumerate(chunks, start=1):
            p = chunks_dir / f"chunk_{i:04d}.list"
            with p.open("w", encoding="utf-8") as fh:
                fh.write("\n".join(chunk) + ("\n" if chunk else ""))
            paths.append(p)
        return paths

    def _process_chunk(tsv_file: Path, pdb_file_path: Path,
                       chunk_list_path: Path, out_fasta_path: Path, out_status_path: Path):
        """Your original per-file loop, restricted to one chunk."""
        uniprot_index = build_uniprot_index(tsv_file)
        current, obsolete = 0, 0

        with chunk_list_path.open("r", encoding="utf-8") as cf, \
             out_fasta_path.open("w", encoding="utf-8") as fout:
            for pdb_file_name in (ln.strip() for ln in cf if ln.strip()):
                if "-trim.pdb" not in pdb_file_name:
                    continue
                pdb_uniprot_id = pdb_file_name.split("AF-")[1].split("-")[0]
                try:
                    e = uniprot_index[pdb_uniprot_id]
                    # Prefer gene unless UC, then fall back to protein (already sanitized upstream)
                    symbol = e.get("gene", "UC")
                    if symbol == "UC":
                        symbol = e.get("protein", "UC")
                    header = f">{e['superkingdom']}|{e['uniprot_id']}|{symbol}|{e['organism']}|{e['rank']}\n"
                    sequence = extract_protein_sequence_with_explicit_gaps(str(pdb_file_path / pdb_file_name))
                    fout.write(header)
                    fout.write(sequence + "\n")
                    current += 1
                except KeyError:
                    obsolete += 1

        with out_status_path.open("w", encoding="utf-8") as sf:
            sf.write(f"Number of total UniProt entries: {current + obsolete}\n")
            sf.write(f"Number of current UniProt entries: {current}\n")
            sf.write(f"Number of obsolete UniProt entries: {obsolete}\n")

    # -------------------- paths & dirs -----------------------------
    SCRIPT_PATH = Path(__file__).resolve()
    LOGS_DIR = OUT_DIR / "logs_trimseq"
    CHUNKS_DIR = OUT_DIR / "chunks_trimseq"
    RESULTS_DIR = OUT_DIR / "trimseq_results"
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Are we inside an array task? (LSF sets LSB_JOBINDEX=1..N)
    job_index = os.environ.get("LSB_JOBINDEX")

    if job_index:
        # -------------------- RUN MODE (compute node) --------------------
        i = int(job_index)
        chunk_file = CHUNKS_DIR / f"chunk_{i:04d}.list"
        out_fasta = RESULTS_DIR / f"trim_blastp_sequences.chunk_{i:04d}.fa"
        out_status = RESULTS_DIR / f"obsolete_status.chunk_{i:04d}.txt"
        if not chunk_file.exists():
            raise FileNotFoundError(f"Chunk list not found: {chunk_file}")
        _process_chunk(TSV_FILE, PDB_FILE_PATH, chunk_file, out_fasta, out_status)
        print(f"[{JOB_NAME}] Finished chunk {i}: wrote {out_fasta.name}, {out_status.name}")

    else:
        # -------------------- SUBMIT MODE (login node) -------------------
        # 1) Gather and filter the PDB file list
        pdb_file_list = sorted([fn for fn in os.listdir(PDB_FILE_PATH) if fn.endswith("-trim.pdb")])
        if not pdb_file_list:
            print("No matching -trim.pdb files found; nothing to do.")
            raise SystemExit(0)

        # 2) Split into N chunks and write lists
        chunks = _chunkify(pdb_file_list, N_CHUNKS)
        chunk_files = _write_chunk_lists(CHUNKS_DIR, chunks)
        n_chunks = len(chunk_files)
        print(f"Prepared {n_chunks} chunk lists in {CHUNKS_DIR}")

        # 3) Submit a single LSF array (each task processes one chunk)
        array_spec = f'"{JOB_NAME}[1-{n_chunks}]"'
        bsub_cmd = (
            f'bsub -q {QUEUE} -J {array_spec} '
            f'-o "{LOGS_DIR}/%J_%I.out" -e "{LOGS_DIR}/%J_%I.err" '
            f'-cwd "{OUT_DIR}" '
            f'python "{SCRIPT_PATH}"'
        )
        print("Submitting LSF array:")
        print(bsub_cmd)
        subprocess.run(bsub_cmd, shell=True, check=True)

        print("\nSubmitted array. When all tasks finish, merge FASTAs with:")
        print(f'  cat "{RESULTS_DIR}"/trim_blastp_sequences.chunk_*.fa > "{OUT_DIR}/trim_blastp_sequences.fa"')
