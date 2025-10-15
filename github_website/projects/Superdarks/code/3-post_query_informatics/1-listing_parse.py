import re
from pathlib import Path
from typing import Dict, Union, Optional

def parse_ls_uniprot_scores(file_path: Union[str, Path]) -> Dict[str, Dict[str, Optional[float]]]:
    """
    Parse an 'ls -l' output file where the final token (filename) encodes:
      rank--coverage-tm--AF-<UNIPROT>-...

    Rules implemented:
      1) strip line and split on spaces
      2) take the last token (filename)
      3) split filename on '--'
      4) tokens:
         4.1) token1  -> rank
         4.2) token2  -> split on '-' : coverage (token1), tm (token2)
         4.3) token3  -> split on 'AF-' then split on '-' : first token is UniProt ID
      5) save {rank, coverage, tm} under key = UniProt ID
      6) return the dictionary

    Returns:
      { 'A0A031HJE0': {'rank': 1, 'coverage': 0.7, 'tm': 0.51}, ... }
    """
    file_path = Path(file_path)
    result: Dict[str, Dict[str, Optional[float]]] = {}

    with file_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("total "):
                continue

            parts = line.split()
            if not parts:
                continue

            last_token = parts[-1]  # filename with embedded scores, per your convention

            # Expect something like:
            # 1--0.7-0.51--AF-A0A031HJE0-F1-model_v4--RHO-...--OK-full.pdb
            segs = last_token.split("--")
            if len(segs) < 3:
                # Not in the expected format; skip
                continue

            # 4.1) rank
            rank_str = segs[0]

            # 4.2) coverage and tm (e.g., "0.7-0.51")
            cov_tm = segs[1].split("-")
            if len(cov_tm) < 2:
                continue
            cov_str, tm_str = cov_tm[0], cov_tm[1]

            # 4.3) UniProt ID from "AF-<ID>-..."
            t3 = segs[2]
            uniprot_id: Optional[str] = None
            if "AF-" in t3:
                after_af = t3.split("AF-", 1)[1]     # e.g., "A0A031HJE0-F1-model_v4"
                uniprot_id = after_af.split("-", 1)[0]
            else:
                # Fallback: try to locate AF-<ID>- pattern anyway
                m = re.search(r"AF-([A-Z0-9]{6,10})-", t3)
                if m:
                    uniprot_id = m.group(1)

            # 4.4) type (OK, NTC, etc...)
            type1 = segs[-1] if segs else ""
            type2 = ""
            if len(segs) > 1:
                s2 = type1.split("-")
                if s2:
                    type2 = s2[0]
            type_class = type2

            if not uniprot_id:
                continue  # cannot identify ID -> skip

            # Convert numbers
            def _to_int_or_float(s: str) -> Optional[float]:
                try:
                    # rank often an integer; keep as int if possible
                    return int(s)
                except Exception:
                    try:
                        return float(s)
                    except Exception:
                        return None

            def _to_float(s: str) -> Optional[float]:
                try:
                    return float(s)
                except Exception:
                    return None

            rank_val = _to_int_or_float(rank_str)
            cov_val = _to_float(cov_str)
            tm_val = _to_float(tm_str)

            result[uniprot_id] = {
                "rank": rank_val,
                "coverage": cov_val,
                "tm": tm_val,
                "type": type_class
            }

    return result
