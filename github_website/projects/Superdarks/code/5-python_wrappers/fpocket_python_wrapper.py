# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.

from __future__ import annotations
from pathlib import Path
import shutil, subprocess, tempfile, os

def _which_fpocket(fpocket_exe: str | None = None) -> str:
    if fpocket_exe:
        return fpocket_exe
    for name in ("fpocket", "fpocket.exe"):
        p = shutil.which(name)
        if p:
            return p
    raise FileNotFoundError("fpocket not found on PATH. Pass fpocket_exe='...' or add it to PATH.")

def _gather_outputs(search_dirs: list[Path], expected_name: str | None = None) -> list[Path]:
    seen: set[Path] = set()
    outs: list[Path] = []
    for root in search_dirs:
        if not root or not root.exists():
            continue
        # Prefer exact expected dir if we know it
        if expected_name:
            cand = root / expected_name
            if cand.is_dir() and cand not in seen:
                seen.add(cand); outs.append(cand)
        # Also collect any *_out dirs (batch/list mode)
        for p in root.glob("*_out"):
            if p.is_dir() and p not in seen:
                seen.add(p); outs.append(p)
    return sorted(outs)

def run_fpocket(
    *,
    pdb_file: str | os.PathLike | None = None,     # -f
    pdb_list_file: str | os.PathLike | None = None,# -F (plain text file with paths)
    pdb_list: list[str | os.PathLike] | None = None, # Python list -> we create tmp list

    # Options (only passed if not None)
    min_radius: float | None = None,   # -m
    max_radius: float | None = None,   # -M
    model_index: int | None = None,    # -l
    cluster_method: str | None = None, # -C
    distance_metric: str | None = None,# -e
    min_alpha_spheres: int | None = None, # -i
    apolar_rule: int | None = None,    # -A
    hier_distance: float | None = None,# -D
    max_apolar_ratio: float | None = None, # -p
    mc_iterations: int | None = None,  # -v
    grid_discretization: int | None = None, # -b
    condensed: bool = False,           # -d
    restrict_ligand: str | None = None,# -r
    pocket_residues: str | None = None,# -P
    amber_top: str | os.PathLike | None = None, # -y
    energy: bool = False,              # -x
    delete_chains: str | None = None,  # -c
    keep_chains: str | None = None,    # -k
    ligand_chain: str | None = None,   # -a
    write_format: str | None = None,   # -w

    # Execution
    workdir: str | os.PathLike | None = None,   # default: PDB parent (if single file)
    fpocket_exe: str | None = None,
    capture_output: bool = True,
    check: bool = True,
) -> dict:
    exe = _which_fpocket(fpocket_exe)

    if sum(x is not None for x in (pdb_file, pdb_list_file, pdb_list)) != 1:
        raise ValueError("Provide exactly one of: pdb_file OR pdb_list_file OR pdb_list.")

    # If a Python list is provided, write a temporary list file for -F
    temp_list_path = None
    if pdb_list is not None:
        tf = tempfile.NamedTemporaryFile(delete=False, suffix=".lst", mode="w", encoding="utf-8")
        for item in pdb_list:
            tf.write(str(item) + "\n")
        tf.close()
        temp_list_path = tf.name
        pdb_list_file = temp_list_path

    # Determine expected output directory name for single-file mode
    expected_out_name = None
    pdb_parent = None
    if pdb_file is not None:
        pdb_path = Path(pdb_file).resolve()
        pdb_parent = pdb_path.parent
        expected_out_name = f"{pdb_path.stem}_out"

    # Default workdir: if not provided and single-file mode, use the PDB parent
    wd = Path(workdir).resolve() if workdir else (pdb_parent if pdb_parent else Path.cwd())
    wd.mkdir(parents=True, exist_ok=True)

    cmd: list[str] = [exe]
    if pdb_file is not None:
        cmd += ["-f", str(Path(pdb_file).resolve())]
    elif pdb_list_file is not None:
        cmd += ["-F", str(Path(pdb_list_file).resolve())]

    def add(flag, val):
        if val is not None:
            cmd.extend([flag, str(val)])
    add("-m", min_radius)
    add("-M", max_radius)
    add("-l", model_index)
    add("-C", cluster_method)
    add("-e", distance_metric)
    add("-i", min_alpha_spheres)
    add("-A", apolar_rule)
    add("-D", hier_distance)
    add("-p", max_apolar_ratio)
    add("-v", mc_iterations)
    add("-b", grid_discretization)
    add("-r", restrict_ligand)
    add("-P", pocket_residues)
    add("-y", amber_top)
    add("-c", delete_chains)
    add("-k", keep_chains)
    add("-a", ligand_chain)
    add("-w", write_format)
    if condensed:
        cmd.append("-d")
    if energy:
        cmd.append("-x")

    try:
        if capture_output:
            proc = subprocess.run(cmd, cwd=str(wd), text=True, capture_output=True, check=False)
        else:
            proc = subprocess.run(cmd, cwd=str(wd), text=True, check=False)
        if check and proc.returncode != 0:
            raise RuntimeError(
                f"fpocket failed (rc={proc.returncode}).\nCommand: {' '.join(cmd)}\nSTDERR:\n{proc.stderr}"
            )
    finally:
        # Clean temp list if we created it
        if temp_list_path:
            try: Path(temp_list_path).unlink()
            except Exception: pass

    # Search in both the run directory and the PDB's parent (if different)
    search_dirs = [wd]
    if pdb_parent and pdb_parent != wd:
        search_dirs.append(pdb_parent)
    output_dirs = _gather_outputs(search_dirs, expected_out_name)

    return {
        "cmd": cmd,
        "workdir": str(wd),
        "expected_out_name": expected_out_name,
        "returncode": proc.returncode,
        "stdout": getattr(proc, "stdout", None),
        "stderr": getattr(proc, "stderr", None),
        "output_dirs": output_dirs,
    }


pdb_file = "./TMEM184C-1--0.86-0.53--AF-Q9NVA4-F1-model_v4--RHO-1F88-1-348-A--OK-edited.pdb"
res = run_fpocket(
    pdb_file=pdb_file,
    min_radius=3.8,   # -m

)
