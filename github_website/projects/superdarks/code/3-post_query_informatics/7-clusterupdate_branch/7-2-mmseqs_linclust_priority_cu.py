# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Priority reps with MMseqs2 (Proteins): Archaea → Bacteria → Eukaryota
Safe DB naming for clusterupdate (no in-place reuse).

Requires: mmseqs on PATH
"""

import subprocess, shutil, sys, tempfile
from pathlib import Path

# -------- CONFIG --------
archaea_faa   = Path("/Users/danisom/Desktop/local_informatics/2023.11.12.trim_blastp_sequences_final.archaea.fa")
bacteria_faa  = Path("/Users/danisom/Desktop/local_informatics/2023.11.12.trim_blastp_sequences_final.bacteria.fa")
euks_faa      = Path("/Users/danisom/Desktop/local_informatics/2023.11.12.trim_blastp_sequences_final.eukaryota.fa")

# Base DB names
A_db,  A_clu        = "A",          "A_clu"
B_db                = "B"
AB_db               = "AB"          # concat(A,B)
ABU_db, ABU_clu     = "AB_updated", "AB_clu"     # outputs from clusterupdate step 1
E_db                = "E"
ABE_db              = "ABE"         # concat(AB_updated, E)
ABEU_db, ABEU_clu   = "ABE_updated","ABE_clu"    # outputs from clusterupdate step 2

# linclust knobs (widely-supported)
MIN_SEQ_ID = 0.50      # --min-seq-id
COV_FRAC   = 0.50      # -c
COV_MODE   = 2         # 1=target, 0=query, 2=both
THREADS    = "14"

REP_DB     = "ABE_repDB"
REP_FASTA  = Path("ABE_reps.fasta")
# ------------------------


def run(cmd: list[str], check=True):
    print(">", " ".join(map(str, cmd)), flush=True)
    return subprocess.run(cmd, check=check)

def need(*files: Path):
    for f in files:
        if not Path(f).exists():
            sys.exit(f"ERROR: Missing input FASTA: {f}")

def rmdb(base: str):
    """Remove an MMseqs DB if present (ignore errors)."""
    if Path(f"{base}.dbtype").exists():
        run(["mmseqs", "rmdb", base], check=False)

def concat_dbs(db1: str, db2: str, outdb: str):
    # (re)create concat outputs (and header DBs)
    rmdb(outdb)
    rmdb(f"{outdb}_h")
    run(["mmseqs", "concatdbs", db1, db2, outdb])
    run(["mmseqs", "concatdbs", f"{db1}_h", f"{db2}_h", f"{outdb}_h"])

def clusterupdate_6(old_db: str, new_superset_db: str, old_clu: str,
                    out_db_updated: str, out_clu_updated: str, tmp: Path):
    # Only remove outputs (never inputs!)
    rmdb(out_db_updated)
    rmdb(out_clu_updated)
    run(["mmseqs", "clusterupdate",
         old_db, new_superset_db, old_clu,
         out_db_updated, out_clu_updated, str(tmp)])

def main():
    if not shutil.which("mmseqs"):
        sys.exit("ERROR: `mmseqs` not found on PATH.")
    need(archaea_faa, bacteria_faa, euks_faa)

    tmp_dir = Path(tempfile.mkdtemp(prefix="mmseqs_tmp_"))

    # 1) createdb (force amino-acid DBs)
    run(["mmseqs", "createdb", "--dbtype", "1", str(archaea_faa),  A_db])
    run(["mmseqs", "createdb", "--dbtype", "1", str(bacteria_faa), B_db])
    run(["mmseqs", "createdb", "--dbtype", "1", str(euks_faa),     E_db])

    # 2) cluster Archaea → A_clu
    rmdb(A_clu)
    run([
        "mmseqs", "linclust",
        A_db, A_clu, str(tmp_dir),
        "--min-seq-id", str(MIN_SEQ_ID),
        "-c", str(COV_FRAC),
        "--cov-mode", str(COV_MODE),
        "--threads", THREADS,
    ])

    # 3) concat A + B => AB ; update A_clu → AB_updated/AB_clu
    concat_dbs(A_db, B_db, AB_db)
    clusterupdate_6(A_db, AB_db, A_clu, ABU_db, ABU_clu, tmp_dir)

    # 4) concat AB_updated + E => ABE ; update AB_clu → ABE_updated/ABE_clu
    concat_dbs(ABU_db, E_db, ABE_db)
    clusterupdate_6(ABU_db, ABE_db, ABU_clu, ABEU_db, ABEU_clu, tmp_dir)

    # 5) reps → FASTA (use the *updated* pair)
    rmdb(REP_DB)
    run(["mmseqs", "createsubdb", ABEU_clu, ABEU_db, REP_DB])
    run(["mmseqs", "convert2fasta", REP_DB, str(REP_FASTA)])

    print("\nDone.")
    print(f"  Representatives DB : {REP_DB}")
    print(f"  Representatives FASTA → {REP_FASTA.resolve()}")

if __name__ == "__main__":
    main()
