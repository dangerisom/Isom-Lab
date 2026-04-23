# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.

from pathlib import Path
from os import sep as _SEP
from collections import OrderedDict
import os
from typing import Optional
from compGeometry import Vertex4D, distance, circumSphere
from pdbFile import PseudoAtom


# -- Example usage ------------------------------------------
# # A) Direct atoms (already CA atoms)
# ca_list = [...]  # list[Atom]
# res = grab_ssco(
#     ca_atoms=ca_list,
#     name="my_model",               # used for output filenames if writing
#     dst_directory="/tmp/out",
#     write_combined=True,
#     return_mode="metrics",         # or "atoms" / "write"
# )

# # B) Full file path
# res = grab_ssco(
#     pdb_file="/data/proteins/FOO123.pdb",  # full path
#     dst_directory="/tmp/out",
#     write_helix=True, write_sheet=True, write_combined=True,
#     return_mode="write",
# )

# # C) Old style (unchanged)
# res = grab_ssco(
#     src_directory="/data/proteins",
#     pdb_file="FOO123.pdb",
#     dst_directory="/tmp/out",
#     write_combined=True,
# )


def _is_atom_like(obj) -> bool:
    """Duck-typing: minimal fields we rely on."""
    return all(hasattr(obj, a) for a in ("x", "y", "z", "residue_sequence_number", "temperature_factor"))

def _compute_ca_certainty(ca_atoms) -> float:
    """Fraction with pLDDT (temperature_factor) > 70, safe for empty."""
    if not ca_atoms:
        return 0.0
    sure = sum(1 for a in ca_atoms if getattr(a, "temperature_factor", 0.0) > 70.0)
    return sure / len(ca_atoms)

def _normalize_file_inputs(src_directory, pdb_file):
    """
    Accept either:
      - (src_directory, pdb_file) as before, or
      - pdb_file as a full path (src_directory may be None)
    Return (src_dir_str_with_sep, pdb_file_name, base_stem)
    """
    if pdb_file is None:
        return None, None, None

    # If user passed a full path in pdb_file, prefer that.
    p = Path(pdb_file)
    if p.exists() and p.is_file():
        src_dir_str = str(p.parent)
        if not src_dir_str.endswith(_SEP):
            src_dir_str += _SEP
        return src_dir_str, p.name, p.stem

    # Fallback to the old style (src_directory + pdb_file name)
    src_dir_str = str(src_directory) if src_directory is not None else ""
    if src_dir_str and not src_dir_str.endswith(_SEP):
        src_dir_str += _SEP
    return src_dir_str, str(pdb_file), Path(pdb_file).stem

def open_pdb(pdb_file_path, pdb_file_name, pdb_format="pdb", zip_status=0):
    # Create an instance of a PDBfile (really protein) object
    if pdb_format == "mmCIF":
        import pdbFile_cif
        pdb_file_object = pdbFile_cif.PDBfile(pdb_file_path, pdb_file_name, zip_status=zip_status)
    else:
        import pdbFile
        pdb_file_object = pdbFile.PDBfile(pdb_file_path, pdb_file_name, zip_status=zip_status)
    return pdb_file_object

def get_ca_atoms(residue_dict, origin=""):
    residue_keys = sorted(residue_dict)
    i, caAtoms, caCertainty = 0, [], 0.0
    while i < len(residue_keys):
        key = residue_keys[i]
        try:
            caAtom = residue_dict[key].atoms["CA"]
            if origin:
                caAtom.symbol = origin
            caAtoms.append(caAtom)
            # AF certainty proxy
            if caAtom.temperature_factor > 70.0:
                caCertainty += 1.0
        except:
            pass
        i += 1
    return (caAtoms, caCertainty/len(caAtoms) if caAtoms else 0.0)

def get_4D_vertices(atom_list):
    v = []
    for atom in atom_list:
        x, y, z = atom.x, atom.y, atom.z
        u = atom.x**2 + atom.y**2 + atom.z**2
        vertex = Vertex4D((x, y, z, u), setC=1)
        vertex.data = atom
        vertex.id = atom.atom_serial
        v.append(vertex)
    return v

# ---------------------------
# Helpers for (3) Optional writes
# ---------------------------
def _residue_key_for_sort(ca):
    """Robustly fetch a sortable key for a CA's residue."""
    # Preferred: ca.residue.key (as requested)
    try:
        key = ca.residue.key
        # If key itself isn't orderable, str() fallback
        _ = (key < key)  # will raise if not comparable; we don't care about result
        return key
    except Exception:
        pass
    # Fallbacks: sequence number, then atom serial
    seq = getattr(ca, "residue_sequence_number", None)
    if seq is not None:
        return (seq, getattr(ca, "atom_serial", 0))
    return getattr(ca, "atom_serial", 0)

def _sorted_ca_dict_by_residue_key(ca_iterable):
    """
    Return OrderedDict sorted by ca.residue.key (fallbacks applied if missing).
    Keys are the residue keys; values are the CA atoms.
    """
    tmp = {}
    for ca in ca_iterable:
        key = getattr(getattr(ca, "residue", object()), "key", None)
        if key is None:
            key = _residue_key_for_sort(ca)
        tmp[key] = ca  # de-duplicate by residue key
    try:
        items = sorted(tmp.items(), key=lambda kv: kv[0])
    except Exception:
        items = sorted(tmp.items(), key=lambda kv: str(kv[0]))
    return OrderedDict(items)


def grab_ssco(
    src_directory=None,
    dst_directory=None,
    pdb_file=None,
    *,
    ca_atoms=None,            
    name: Optional[str] = None,
    pdb_format: str = "pdb",
    zip_status: int = 0,
    write_helix: bool = False,
    write_sheet: bool = False,
    write_combined: bool = False,
    return_mode: str = "metrics",   # "atoms" | "write" | "metrics"
):
    """
    Compute SSCO (structured C-alpha fraction) for one structure and optionally
    write PDBs containing only helical, only sheet(-like), or combined residues.

    INPUT MODES
    -----------
    A) Direct atoms:
        - pass ca_atoms=[Atom, Atom, ...]
        - optional name="basename" for output files
    B) Full file path:
        - pass pdb_file="/path/to/file.pdb", src_directory can be None
    C) Old style:
        - pass src_directory="/path/to/", pdb_file="file.pdb"

    return_mode:
      - "atoms"   -> returns OrderedDict sorted by residue key (combined CA atoms)
      - "write"   -> returns dict describing written outputs (paths + counts)
      - "metrics" -> returns (ssco_score, ca_certainty, n_ca)
    """
    # ─────────────────────────────────────────────
    # Resolve input mode -> (ca_atoms, ca_certainty, n_ca, base_name)
    # ─────────────────────────────────────────────
    base_name = None

    if ca_atoms is not None:
        # Accept a single list of Atom objects (assumed to be CA atoms)
        if not isinstance(ca_atoms, (list, tuple)) or not all(_is_atom_like(a) for a in ca_atoms):
            raise ValueError("ca_atoms must be a list of Atom-like objects with coordinates and residue info.")
        # If caller accidentally passed all atoms, keep only CA
        ca_atoms = [a for a in ca_atoms if getattr(a, "atom_name", "").strip().upper() == "CA" or a.atom_name == " CA " or a.atom_name == "CA"]
        ca_certainty = _compute_ca_certainty(ca_atoms)
        n_ca = len(ca_atoms)
        base_name = name or "ca_input"
    else:
        # File-based modes
        src_dir_str, pdb_file_name, base_name = _normalize_file_inputs(src_directory, pdb_file)
        if not pdb_file_name:
            raise ValueError("Must provide either ca_atoms or a pdb_file (full path or with src_directory).")
        query = open_pdb(src_dir_str, pdb_file_name, pdb_format=pdb_format, zip_status=zip_status)
        ca_atoms, ca_certainty = get_ca_atoms(query.residues)
        n_ca = len(ca_atoms)

    # If no CA atoms, honor return_mode and exit early
    if n_ca == 0:
        if return_mode == "atoms":
            return OrderedDict()
        if return_mode == "write":
            return {
                "helix": {"written": False, "path": None, "n": 0},
                "sheet": {"written": False, "path": None, "n": 0},
                "combined": {"written": False, "path": None, "n": 0},
            }
    # Build vertices
    ca_vertices = get_4D_vertices(ca_atoms)

    # --------------------------
    # 1) Identify helical CA's
    # --------------------------
    helix_map = {}  # seq_no -> CA atom
    for i in range(0, max(0, len(ca_vertices) - 3)):
        v1, v2, v3, v4 = ca_vertices[i], ca_vertices[i+1], ca_vertices[i+2], ca_vertices[i+3]
        psa = PseudoAtom()
        cs_radius = circumSphere(v1, v2, v3, v4, psa, returnSquareDistance=0)[0]
        if 2.75 < cs_radius < 3.25:
            helix_map[v1.data.residue_sequence_number] = v1.data
            helix_map[v2.data.residue_sequence_number] = v2.data
            helix_map[v3.data.residue_sequence_number] = v3.data
            helix_map[v4.data.residue_sequence_number] = v4.data

    # group into runs (gap ≤2)
    helix_runs, run = [], {}
    h_keys = sorted(helix_map)
    for i in range(0, max(0, len(h_keys) - 1)):
        ca1 = helix_map[h_keys[i]]
        ca2 = helix_map[h_keys[i+1]]
        if abs(ca1.residue_sequence_number - ca2.residue_sequence_number) <= 2:
            run[ca1] = None
            run[ca2] = None
        else:
            if run:
                helix_runs.append(list(run))
            run = {}
    if run:
        helix_runs.append(list(run))

    # keep long helix runs
    h_keep_cas = []
    for r in helix_runs:
        if len(r) > 8:
            h_keep_cas.extend(r)
    h_keep_set = set(h_keep_cas)

    # non-helical vertices for sheet detection
    unclassified = [v for v in ca_vertices if v.data not in h_keep_set]

    # --------------------------
    # 2) Identify sheet(-like)
    # --------------------------
    sheet_map = {}  # seq_no -> CA atom
    for i in range(2, max(2, len(unclassified) - 3)):
        vm2 = unclassified[i-2]
        vm1 = unclassified[i-1]
        v1  = unclassified[i]
        vp1 = unclassified[i+1]
        vp2 = unclassified[i+2]
        vp3 = unclassified[i+3]

        sheet_distance = distance(v1, vp3)

        # Count contacts for v1 (skip immediate neighbors)
        contacts = 0
        for v2 in unclassified:
            if v2 in (vm2, vm1, v1, vp1, vp2):
                continue
            if distance(v1, v2) < 8.0:
                contacts += 1

        # Modified 2025.08.24 after bug in logic found
        # I changed contacts and sheet_distance to replicate the original rhod ssco pattern
        """
        1) A bug in ssco.py inflates sheet contacts
        In ssco.py the contact-count loop meant to skip immediate neighbors uses the wrong variable, 
        so neighbors aren’t actually skipped. That inflates the contact count and changes which residues get called “sheet/sheet-like”.
        
        Current (buggy) logic in ssco.py:
        for v2 in query_ca_vertices_unclassified:
            if v1 not in [vm1, vm2, vp1, vp2]:   # <-- should be v2, not v1
                d = distance(v1, v2)
                if d < 8:
                    contacts += 1
        Fix:
        for v2 in query_ca_vertices_unclassified:
            if v2 in (vm2, vm1, v1, vp1, vp2):   # skip immediate neighbors (including self)
                continue
            if distance(v1, v2) < 8.0:
                contacts += 1

        """
        # if contacts > 4 and sheet_distance > 5.0:
        #     sheet_map[v1.data.residue_sequence_number] = v1.data
        if contacts > 0 and sheet_distance > 6.0:
            sheet_map[v1.data.residue_sequence_number] = v1.data

    # group sheet runs (gap ≤1)
    sheet_runs, run = [], {}
    s_keys = sorted(sheet_map)
    for i in range(0, max(0, len(s_keys) - 1)):
        ca1 = sheet_map[s_keys[i]]
        ca2 = sheet_map[s_keys[i+1]]
        if abs(ca1.residue_sequence_number - ca2.residue_sequence_number) <= 1:
            run[ca1] = None
            run[ca2] = None
        else:
            if run:
                sheet_runs.append(list(run))
            run = {}
    if run:
        sheet_runs.append(list(run))

    # drop short runs (≤4)
    sheet_runs = [r for r in sheet_runs if len(r) > 4]

    # keep only runs with inter-run contacts (<5 Å to another run)
    s_keep_runs = []
    for r1 in sheet_runs:
        has_contact = False
        for r2 in sheet_runs:
            if r1 is r2:
                continue
            for ca1 in r1:
                if any(distance(ca1, ca2) < 5.0 for ca2 in r2):
                    has_contact = True
                    break
            if has_contact:
                break
        if has_contact:
            s_keep_runs.append(r1)

    # within-run flags: any longer-range (|Δi|>4) contact < 6 Å?
    s_flags_per_run = []
    for r1 in s_keep_runs:
        flags = {ca: 0 for ca in r1}
        for r2 in s_keep_runs:
            for ca1 in r1:
                for ca2 in r2:
                    if abs(ca1.residue_sequence_number - ca2.residue_sequence_number) > 4:
                        if distance(ca1, ca2) < 6.0:
                            flags[ca1] = 1
        s_flags_per_run.append(flags)

    # trim loose ends (leading/trailing zeros) keeping run order
    trim_nums = []
    for flags in s_flags_per_run:
        cas = list(flags)
        i = 0
        while i < len(cas) and flags[cas[i]] == 0:
            trim_nums.append(cas[i].residue_sequence_number)
            i += 1
        i = len(cas) - 1
        while i >= 0 and flags[cas[i]] == 0:
            trim_nums.append(cas[i].residue_sequence_number)
            i -= 1

    # final keep set for sheets
    s_keep_cas = []
    for r in s_keep_runs:
        if len(r) >= 4:
            for ca in r:
                if ca.residue_sequence_number not in trim_nums:
                    s_keep_cas.append(ca)

    # --------------------------
    # 3) Optional writes (+ return options; combined sorted by residue.key)
    # --------------------------
    dst = Path(dst_directory or ".")
    base = base_name or (Path(pdb_file).stem if pdb_file else "ssco")

    def _write_residues(res_list, out_path: Path):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8") as f:
            for ca in res_list:
                f.write(str(ca.residue))

    write_info = {
        "helix": {"written": False, "path": None, "n": 0},
        "sheet": {"written": False, "path": None, "n": 0},
        "combined": {"written": False, "path": None, "n": 0},
    }

    # Write helix
    if write_helix and h_keep_cas:
        hp = dst / f"{base}-h.pdb"
        _write_residues(h_keep_cas, hp)
        write_info["helix"] = {"written": True, "path": hp, "n": len(h_keep_cas)}

    # Write sheet
    if write_sheet and s_keep_cas:
        sp = dst / f"{base}-s.pdb"
        _write_residues(s_keep_cas, sp)
        write_info["sheet"] = {"written": True, "path": sp, "n": len(s_keep_cas)}

    # Combined (sorted by ca.residue.key)
    combined_sorted_od = _sorted_ca_dict_by_residue_key(list(h_keep_cas) + list(s_keep_cas))
    if write_combined and combined_sorted_od:
        cp = dst / f"{base}-ssco.pdb"
        cp.parent.mkdir(parents=True, exist_ok=True)
        with cp.open("w", encoding="utf-8") as f:
            for _, ca in combined_sorted_od.items():  # sorted order
                f.write(str(ca.residue))
        write_info["combined"] = {"written": True, "path": cp, "n": len(combined_sorted_od)}

    # --------------------------
    # 4) Score & return
    # --------------------------
    total_structured = float(len(h_keep_cas) + len(s_keep_cas))
    ssco_score = (total_structured / n_ca) if n_ca else 0.0

    if return_mode == "atoms":
        return combined_sorted_od
    if return_mode == "write":
        return write_info
    return (ssco_score, ca_certainty, n_ca)