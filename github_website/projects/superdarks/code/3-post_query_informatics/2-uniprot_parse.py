# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.


import re
from typing import Dict, List, Optional

def _collect_lines(entry_text: str, prefix: str) -> List[str]:
    lines = []
    for ln in entry_text.splitlines():
        if ln.startswith(prefix):
            lines.append(ln[len(prefix):].rstrip())
    return lines

def _parse_accessions(entry_text: str) -> Dict[str, List[str]]:
    """Return {"primary": <str or "">, "all": [..]}."""
    ac_lines = _collect_lines(entry_text, "AC   ")
    accs: List[str] = []
    for part in ac_lines:
        # Tokens are like 'P12345; Q9ABC1;'
        for tok in part.split(";"):
            tok = tok.strip()
            if not tok:
                continue
            # Keep only accession-like tokens
            # NEW (6-char OR 10-char A0A…)
            if re.match(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|A0[A-Z0-9]{8})$", tok):
                accs.append(tok)
    primary = accs[0] if accs else ""
    return {"primary": primary, "all": accs}

def _parse_status(entry_text: str) -> Optional[str]:
    # ID   PROT_HUMAN Reviewed; 350 AA.
    m = re.search(r"^ID\s+.+?\s+(Reviewed|Unreviewed);", entry_text, flags=re.M)
    return m.group(1) if m else None

def _parse_id_line(entry_text: str) -> Optional[str]:
    m = re.search(r"^ID\s+.*$", entry_text, flags=re.M)
    return m.group(0) if m else None

def _parse_recommended_protein(entry_text: str) -> Optional[str]:
    # Look for 'DE   RecName: Full=...;' else 'DE   SubName: Full=...;'
    de_text = "\n".join(_collect_lines(entry_text, "DE   "))
    if not de_text:
        return None
    m = re.search(r"RecName:\s*Full=([^;]+);", de_text)
    if m:
        return m.group(1).strip()
    m = re.search(r"SubName:\s*Full=([^;]+);", de_text)
    if m:
        return m.group(1).strip()
    return None

def _parse_primary_gene(entry_text: str) -> Optional[str]:
    # GN   Name=TP53; Synonyms=P53; ORFNames=...;
    gn_text = " ".join(_collect_lines(entry_text, "GN   "))
    if not gn_text:
        return None
    m = re.search(r"Name=([^;]+);", gn_text)
    if not m:
        return None
    val = m.group(1).strip()
    # If multiple names separated by commas/spaces, take the first token
    # (Usually primary is a single symbol anyway.)
    for sep in [",", ";", " "]:
        if sep in val:
            return val.split(sep)[0].strip()
    return val

def _parse_organism(entry_text: str) -> Optional[str]:
    # OS   Homo sapiens (Human).
    os_text = " ".join(_collect_lines(entry_text, "OS   ")).strip()
    if not os_text:
        return None
    # Drop trailing period if present
    return os_text[:-1] if os_text.endswith(".") else os_text

def _parse_lineage(entry_text: str) -> List[str]:
    # OC   Eukaryota; Metazoa; Chordata; ...
    oc_text = " ".join(_collect_lines(entry_text, "OC   "))
    parts = [p.strip() for p in oc_text.split(";") if p.strip()]
    return parts

def _parse_interpro_ids(entry_text: str) -> List[str]:
    # DR   InterPro; IPR000001; Some name.
    ids: List[str] = []
    for ln in entry_text.splitlines():
        if ln.startswith("DR   InterPro;"):
            m = re.search(r"InterPro;\s*(IPR[0-9A-Za-z]+);", ln)
            if m:
                ids.append(m.group(1))
    # deduplicate preserving order
    seen = set()
    uniq = []
    for x in ids:
        if x not in seen:
            seen.add(x)
            uniq.append(x)
    return uniq

def _parse_protein_existence(entry_text: str) -> Optional[str]:
    # PE   1: Evidence at protein level
    m = re.search(r"^PE\s+(\d)\s*:", entry_text, flags=re.M)
    return m.group(1) if m else None

def _parse_sequence(entry_text: str) -> Optional[str]:
    # SQ   SEQUENCE   393 AA;  43962 MW;  B3F3D6C86F018D44 CRC64;
    #      MEEPQSDPSV EPPLSQETF... (letters and spaces, sometimes numbers at EOL)
    sq_start = re.search(r"^SQ\s", entry_text, flags=re.M)
    if not sq_start:
        return None
    # sequence lines follow SQ until '//' line
    seq_lines: List[str] = []
    started = False
    for ln in entry_text.splitlines():
        if ln.startswith("SQ "):
            started = True
            continue
        if not started:
            continue
        if ln.startswith("//"):
            break
        # strip spaces and digits
        seq_lines.append(re.sub(r"[^A-Za-z]", "", ln))
    seq = "".join(seq_lines).upper()
    return seq or None

def parse_uniprot_entry(entry_text: str) -> Dict:
    """
    Parse a single UniProtKB flat-file *entry* (text block ending with '//').
    Returns a mapping {primary_accession: {fields}}.
    Only the recommended protein name and primary gene name are included as 'protein' and 'gene'.
    Other useful fields are added for convenience.
    """
    # 1) Accessions
    acc = _parse_accessions(entry_text)
    primary_accession = acc["primary"]
    all_accessions = acc["all"]

    # 2) Minimal checks
    if not primary_accession:
        return {}

    # 3) Core fields
    status = _parse_status(entry_text)
    id_line = _parse_id_line(entry_text)
    protein_recommended = _parse_recommended_protein(entry_text)
    gene_primary = _parse_primary_gene(entry_text)
    organism = _parse_organism(entry_text)
    lineage_list = _parse_lineage(entry_text)
    superkingdom = lineage_list[0] if lineage_list else None
    interpro_ids = _parse_interpro_ids(entry_text)
    pe_level = _parse_protein_existence(entry_text)
    sequence = _parse_sequence(entry_text)

    # 4) Normalize for TSV friendliness
    lineage = "; ".join(lineage_list) if lineage_list else None
    interpro = ";".join(interpro_ids) if interpro_ids else None

    # 5) Assemble output dictionary
    out: Dict[str, dict] = {}
    out[primary_accession] = {
        "status": status,
        "primary_accession": primary_accession,
        "protein": protein_recommended,   # ONLY recommended name
        "gene": gene_primary,             # ONLY primary gene name
        "organism": organism,
        "lineage": lineage,               # semicolon-separated
        "superkingdom": superkingdom,
        "interpro_ids": interpro,         # semicolon-separated
        "protein_existence": pe_level,
        "sequence_aa": sequence,
        "id_line": id_line,
        # "secondary_accessions": ";".join(all_accessions[1:]) if len(all_accessions) > 1 else None,
    }
    return out
