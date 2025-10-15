#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path

def is_index_present_prefix(prefix: Path) -> bool:
    prefix = Path(prefix)
    return (prefix.parent / (prefix.name + ".dbtype")).exists()


# ---- Foldseek createdb file filters (avoid macOS hidden files, non-structures) ----
FOLDSEEK_FILE_INCLUDE = r".*\.(cif|mmcif|pdb)(\.gz)?$"
FOLDSEEK_FILE_EXCLUDE = r"(^|/)\.|\.DS_Store$"


import os
import sys
import queue
import threading
import subprocess
from shlex import quote
from pathlib import Path
from datetime import datetime
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import shutil  # needed by which_or_die
def is_index_fresh(prefix: Path) -> bool:
    """Return True if <prefix>.idx and .idx.index exist and look up-to-date.
    Criteria: both files exist, .idx > ~100 MB, and .idx mtime >= .dbtype mtime.
    """
    try:
        prefix = Path(prefix)
        dbtype = prefix.with_suffix('.dbtype')
        idx = prefix.with_suffix('.idx')
        idx_index = Path(str(prefix) + '.idx.index')
        if not (dbtype.exists() and idx.exists() and idx_index.exists()):
            return False
        if idx.stat().st_size < 100_000_000:
            return False
        return idx.stat().st_mtime >= dbtype.stat().st_mtime
    except Exception:
        return False

# -------------------- Human-readable options --------------------

ALIGNMENT_OPTIONS = [("0", "local(3Di)"),
                     ("1", "TMalign"),
                     ("2", "local(3Di+AA)")]

COVMODE_OPTIONS   = [("0", "query&target"),
                     ("1", "query"),
                     ("2", "target")]

def _opts_to_strings(opts):
    return [f"{k}: {v}" for k, v in opts]

def _code_from_choice(choice: str, default: int = 0) -> int:
    try:
        return int(choice.split(":", 1)[0].strip())
    except Exception:
        return default

# -------------------- Utilities --------------------

def which_or_die(exe: str) -> str:
    path = shutil.which(exe)
    if path is None:
        raise FileNotFoundError(
            f"'{exe}' not found on PATH. Install Foldseek (e.g., `brew install foldseek`)."
        )
    return path

def quote(s: str) -> str:
    return f'"{s}"' if (" " in s or "\t" in s) else s

from pathlib import Path  # already present above

# --- Database prefix utilities ---
def normalize_db_prefix(p: Path) -> Path:
    """
    Return the Foldseek DB *prefix* (path without .dbtype).
    Accepts a .dbtype file, a prefix, or a folder containing exactly one .dbtype.
    """
    p = Path(p)
    # Case 1: user picked the .dbtype file itself
    if p.is_file() and p.suffix == ".dbtype":
        return p.with_suffix("")
    # Case 2: user picked a directory; allow if exactly one .dbtype inside
    if p.is_dir():
        dbtypes = list(p.glob("*.dbtype"))
        if len(dbtypes) == 1:
            return dbtypes[0].with_suffix("")
        raise ValueError("Please select the .dbtype file (or a folder with exactly one .dbtype).")
    # Case 3: user entered a prefix; accept if <prefix>.dbtype exists
    if (p.parent / (p.name + ".dbtype")).exists():
        return p
    raise ValueError("Invalid database selection: need a .dbtype file or a valid DB prefix.")
    

def is_index_present(db_path: Path) -> bool:
    """True if the selected DB (file/dir/prefix) resolves to a valid prefix with a .dbtype."""
    try:
        pref = normalize_db_prefix(db_path)
        return (pref.parent / (pref.name + ".dbtype")).exists()
    except Exception:
        return False


def is_index_present(db_path: Path) -> bool:
    # Heuristic: typical index-related files in a Foldseek DB folder
    if not db_path.exists() or not db_path.is_dir():
        return False
    patterns = ["*.dbtype", "*.idx", "*.lookup"]
    for pat in patterns:
        if list(db_path.glob(pat)):
            return True
    return False

def find_query_files(q: Path) -> list[Path]:
    if q.is_file():
        return [q]
    if q.is_dir():
        files = [p for p in q.rglob("*") if p.suffix.lower() in {".pdb", ".cif"}]
        return sorted(files)
    return []

def write_command_log(out_dir: Path, command: list[str]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log = out_dir / f"foldseek_run_{ts}.cmdlog.txt"
    with open(log, "w", encoding="utf-8") as f:
        f.write(f"# Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# CWD: {Path.cwd()}\n")
        f.write("# Command:\n")
        f.write(" ".join(quote(x) for x in command) + "\n")

def write_named_command_log(log_path: Path, command: list[str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(f"# Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# CWD: {Path.cwd()}\n")
        f.write("# Command:\n")
        f.write(" ".join(quote(x) for x in command) + "\n")

def _fmt_num_for_fname(x: float | int | str) -> str:
    try:
        s = f"{float(x):g}"
    except Exception:
        s = str(x)
    return s.replace('.', 'p').replace('-', 'm')
def build_search_command(
    query_file: Path,
    db_path: Path,
    out_file: Path,
    tmp_dir: Path,
    threads: int,
    sensitivity: float,
    num_iterations: int | None,
    exhaustive_flag: bool,
    max_seqs: int,
    evalue: float,
    alignment_type: int,
    cov: float,
    cov_mode: int,
    gpu_enabled: bool,
    cluster_search: bool = False,
) -> list[str]:
    cmd = [
        "foldseek", "easy-search",
        str(query_file),
        str(db_path),
        str(out_file),
        str(tmp_dir),
        "-s", str(sensitivity),
        "--max-seqs", str(max_seqs),
        "-e", str(evalue),
        "--alignment-type", str(alignment_type),
        "-c", str(cov),
        "--cov-mode", str(cov_mode),
        "--threads", str(threads),
        "--db-load-mode", "2",
    ]
    if cluster_search:
        cmd.extend(["--cluster-search", "1"])
    if num_iterations is not None:
        cmd.extend(["--num-iterations", str(num_iterations)])
    if exhaustive_flag:
        cmd.append("--exhaustive-search")
    # GPU note: CUDA only; not for Apple Silicon/Metal. We keep it off on macOS by default.
    if gpu_enabled:
        cmd.extend(["--gpu", "1"])
    else:
        cmd.extend(["--gpu", "0"])
    return cmd


def run_cmd_streamed(cmd: list[str], logq: queue.Queue):
    env = os.environ.copy()
    base = env.get('CREATEDB_PAR', '').strip()
    filters = f"--file-include '{FOLDSEEK_FILE_INCLUDE}' --file-exclude '{FOLDSEEK_FILE_EXCLUDE}'" if 'FOLDSEEK_FILE_INCLUDE' in globals() else ""
    env['CREATEDB_PAR'] = (base + ' ' + filters).strip()

    logq.put('[CMD] ' + ' '.join(quote(x) for x in cmd) + '\n')
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
        env=env,
    )
    try:
        for line in proc.stdout:
            logq.put(line)
    finally:
        pass
    rc = proc.wait()
    if rc == 0:
        logq.put('[INFO] done.')
    else:
        logq.put(f'[ERROR] Exit code {rc}')

class FoldseekGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Foldseek Runner (macOS)")
        self.geometry("980x660")

        # state
        self.worker: threading.Thread | None = None
        self.log_queue: queue.Queue = queue.Queue()

        # Defaults (Apple Silicon-friendly)
        self.query_path = tk.StringVar(value=str(Path.home()))
        self.db_path = tk.StringVar(value=str(Path.home()))
        self.output_dir = tk.StringVar(value="/Users/danisom/Desktop/foldseek_results/")
        self.tmp_dir = tk.StringVar(value="")  # auto default under DB
        self.sensitivity = tk.DoubleVar(value=9.5)
        self.num_iterations = tk.StringVar(value="1")
        self.max_seqs = tk.IntVar(value=1000)
        self.evalue = tk.DoubleVar(value=1.0)
        # underlying ints still exist
        self.alignment_type = tk.IntVar(value=0)  # 0 local(3Di), 1 TMalign, 2 local(3Di+AA)
        self.cov_mode = tk.IntVar(value=0)
        # display strings for the comboboxes
        align_default = next((f"{k}: {v}" for k, v in ALIGNMENT_OPTIONS
                              if int(k) == self.alignment_type.get()), "1: TMalign")
        covmode_default = next((f"{k}: {v}" for k, v in COVMODE_OPTIONS
                                if int(k) == self.cov_mode.get()), "0: query&target")
        self.alignment_choice = tk.StringVar(value=align_default)
        self.cov_mode_choice  = tk.StringVar(value=covmode_default)

        self.coverage = tk.DoubleVar(value=0.0)
        self.threads = tk.IntVar(value=max(os.cpu_count() - 2, 2))
        self.gpu_enabled = tk.BooleanVar(value=False)  # CUDA-only; OFF on macOS
        self.cluster_search = tk.BooleanVar(value=False)  # --cluster-search off by default
        self.ensure_index_flag = tk.BooleanVar(value=True)
        self.exhaustive_flag = tk.BooleanVar(value=False)

        self._build_ui()
        self.after(100, self._poll_log_queue)

    def _build_ui(self):
        pad = {"padx": 8, "pady": 6}

        # Paths frame
        frm_paths = ttk.LabelFrame(self, text="Paths")
        frm_paths.pack(fill="x", **pad)

        self._row_path_picker(frm_paths, "Query file or folder:", self.query_path, self._pick_query, 0)
        self._row_path_picker(frm_paths, "Database prefix (.dbtype or prefix):", self.db_path, self._pick_db, 1)
        self._row_path_picker(frm_paths, "Output folder:", self.output_dir, self._pick_out, 2)
        self._row_path_picker(frm_paths, "Temp folder (optional):", self.tmp_dir, self._pick_tmp, 3)

        # Options frame
        frm_opts = ttk.LabelFrame(self, text="Search Options")
        frm_opts.pack(fill="x", expand=True, **pad)

        # Columns:
        # 0: label A    (fixed)
        # 1: control A  (expands)
        # 2: label B    (fixed)
        # 3: control B  (expands)
        frm_opts.grid_columnconfigure(0, weight=0)
        frm_opts.grid_columnconfigure(1, weight=1, minsize=160)  # first input column
        frm_opts.grid_columnconfigure(2, weight=0)
        frm_opts.grid_columnconfigure(3, weight=1, minsize=160)  # second input column


        self._row_entry(frm_opts, "Iterations (--num-iterations, empty=None):",
                        self.num_iterations, 1, 2, width=8)
        self._row_entry(frm_opts, "Sensitivity (-s):", self.sensitivity, 1, 0, width=8)
        # (removed) redundant --exhaustive-search checkbox
        self._row_entry(frm_opts, "Max seqs (--max-seqs):", self.max_seqs, 2, 0, width=8)
        self._row_entry(frm_opts, "E-value (-e):", self.evalue, 2, 2, width=8)

        # ----- Row 3: Alignment type (dropdown with labels) + Coverage -----
        ttk.Label(frm_opts, text="Alignment type:")\
            .grid(row=3, column=0, sticky="w", padx=8, pady=6)

        cb_align = ttk.Combobox(
            frm_opts,
            textvariable=self.alignment_choice,
            values=_opts_to_strings(ALIGNMENT_OPTIONS),
            state="readonly", width=22
        )
        cb_align.grid(row=3, column=1, sticky="ew", padx=4, pady=6)

        ttk.Label(frm_opts, text="Coverage (-c):")\
            .grid(row=3, column=2, sticky="w", padx=8, pady=6)
        ttk.Entry(frm_opts, textvariable=self.coverage, width=8)\
            .grid(row=3, column=3, sticky="ew", padx=4, pady=6)

        # ----- Row 4: Cov mode (dropdown with labels) + Threads -----
        ttk.Label(frm_opts, text="Cov mode (--cov-mode):")\
            .grid(row=4, column=0, sticky="w", padx=8, pady=6)

        cb_covmode = ttk.Combobox(
            frm_opts,
            textvariable=self.cov_mode_choice,
            values=_opts_to_strings(COVMODE_OPTIONS),
            state="readonly", width=22
        )
        cb_covmode.grid(row=4, column=1, sticky="ew", padx=4, pady=6)

        ttk.Label(frm_opts, text="Threads (--threads):")\
            .grid(row=4, column=2, sticky="w", padx=8, pady=6)
        ttk.Entry(frm_opts, textvariable=self.threads, width=8)\
            .grid(row=4, column=3, sticky="ew", padx=4, pady=6)

        self._row_check(frm_opts, "Enable GPU (CUDA only; not for macOS Metal)",
                        self.gpu_enabled, 5, 0)
        self._row_check(frm_opts, "Ensure DB index (createindex if missing)", self.ensure_index_flag, 5, 2)
        self._row_check(frm_opts, "Exhaustive search (--exhaustive-search)", self.exhaustive_flag, 6, 2)
        self._row_check(frm_opts, "Report all cluster members (--cluster-search)", self.cluster_search, 6, 0)

        # Buttons
        frm_btn = ttk.Frame(self)
        frm_btn.pack(fill="x", **pad)
        ttk.Button(frm_btn, text="Create Index Now", command=self._create_index).pack(side="left", padx=4)
        ttk.Button(frm_btn, text="Run Search", command=self._run_search).pack(side="left", padx=4)
        ttk.Button(frm_btn, text="Open Output Folder", command=self._open_output).pack(side="left", padx=4)

        # Log area
        frm_log = ttk.LabelFrame(self, text="Log")
        frm_log.pack(fill="both", expand=True, **pad)
        self.txt = tk.Text(frm_log, height=16, wrap="word")
        self.txt.pack(fill="both", expand=True)
        self.txt.insert("end", "Ready.\n")

    # ---------- UI helpers ----------

    def _row_path_picker(self, parent, label, var, picker, row):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=8, pady=6)
        entry = ttk.Entry(parent, textvariable=var, width=80)
        entry.grid(row=row, column=1, sticky="we", padx=4)
        ttk.Button(parent, text="Browse…", command=picker).grid(row=row, column=2, padx=4)
        parent.grid_columnconfigure(1, weight=1)

    def _row_entry(self, parent, label, var, row, col, width=10):
        ttk.Label(parent, text=label).grid(row=row, column=col, sticky="w", padx=8, pady=6)
        e = ttk.Entry(parent, textvariable=var, width=width)
        e.grid(row=row, column=col+1, sticky="ew", padx=4, pady=6)  # expands with column 1 or 3

    def _row_dropdown(self, parent, label, var, values, row, col, tooltip=None):
        ttk.Label(parent, text=label).grid(row=row, column=col, sticky="w", padx=8, pady=6)
        combo = ttk.Combobox(parent, textvariable=var, values=values, state="readonly")
        combo.grid(row=row, column=col+1, sticky="ew", padx=4, pady=6)

    def _row_check(self, parent, label, var, row, col):
        ttk.Checkbutton(parent, text=label, variable=var)\
            .grid(row=row, column=col, columnspan=2, sticky="w", padx=8, pady=6)

    # ---------- Path pickers ----------

    def _pick_query(self):
        # Let user choose either file or folder
        resp = messagebox.askyesno("Choose", "Select YES for a FILE, NO for a FOLDER.")
        if resp:
            p = filedialog.askopenfilename(title="Choose a query file (.pdb/.cif)",
                                           filetypes=[("Structure", "*.pdb *.cif"), ("All", "*.*")])
            if p:
                self.query_path.set(p)
        else:
            d = filedialog.askdirectory(title="Choose a folder containing .pdb/.cif")
            if d:
                self.query_path.set(d)

    def _pick_db(self):
        p = filedialog.askopenfilename(title="Choose Foldseek database .dbtype",
                                       filetypes=[("Foldseek DB index", "*.dbtype"), ("All files", "*.*")])
        if p:
            self.db_path.set(p)
            # auto set tmp dir inside DB prefix parent if empty
            try:
                pref = normalize_db_prefix(Path(p))
                if not self.tmp_dir.get().strip():
                    self.tmp_dir.set(str(pref.parent / "tmp"))
            except Exception:
                pass

    def _pick_out(self):
        d = filedialog.askdirectory(title="Choose output folder")
        if d:
            self.output_dir.set(d)

    def _pick_tmp(self):
        d = filedialog.askdirectory(title="Choose temp folder")
        if d:
            self.tmp_dir.set(d)

    # ---------- Actions ----------

    def _create_index(self):
        try:
            db = normalize_db_prefix(Path(self.db_path.get()))
            tmp = Path(self.tmp_dir.get() or (db / "tmp"))
            tmp.mkdir(parents=True, exist_ok=True)

            if is_index_present(db):
                self._log("[INFO] Database index appears to exist. Skipping.\n")
                return

            cmd = ["foldseek", "createindex", str(db), str(tmp), "--threads", str(self.threads.get())]
            self._run_background([cmd], note="Creating index…")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _open_output(self):
        out = Path(self.output_dir.get())
        out.mkdir(parents=True, exist_ok=True)
        if sys.platform == "darwin":
            subprocess.Popen(["open", str(out)])
        else:
            subprocess.Popen(["xdg-open", str(out)])

    def _run_search(self):
        try:
            qpath = Path(self.query_path.get())
            db = normalize_db_prefix(Path(self.db_path.get()))
            outdir = Path(self.output_dir.get())
            tmp = Path(self.tmp_dir.get() or (db / "tmp"))

            if not qpath.exists():
                raise FileNotFoundError("Query path does not exist.")
            if not is_index_present_prefix(db):
                raise FileNotFoundError("Database index (.dbtype) not found for the selected prefix.")
            outdir.mkdir(parents=True, exist_ok=True)
            tmp.mkdir(parents=True, exist_ok=True)

            # Ensure index if requested
            cmds = []
            if self.ensure_index_flag.get() and not is_index_present(db):
                cmds.append(["foldseek", "createindex", str(db), str(tmp), "--threads", str(self.threads.get())])

            # Resolve queries
            files = find_query_files(qpath)
            if not files:
                raise FileNotFoundError("No .pdb/.cif files found at the query path.")

            # Convert visible dropdown choices -> integer codes
            self.alignment_type.set(_code_from_choice(self.alignment_choice.get(), default=0))
            self.cov_mode.set(_code_from_choice(self.cov_mode_choice.get(), default=0))

            # Build commands for each file
            for f in files:
                base = f.stem.replace(" ", "_")
                out_file = outdir / f"{base}_{'exhaustive' if bool(self.exhaustive_flag.get()) else 'easy'}_s{_fmt_num_for_fname(self.sensitivity.get())}_e{_fmt_num_for_fname(self.evalue.get())}.m8"

                # Interpret iterations (empty -> None)
                iters = self.num_iterations.get().strip()
                iters_val = int(iters) if iters != "" else None

                cmd = build_search_command(
                    query_file=f,
                    db_path=db,
                    out_file=out_file,
                    tmp_dir=tmp,
                    threads=int(self.threads.get()),
                    sensitivity=float(self.sensitivity.get()),
                    num_iterations=iters_val,
                    exhaustive_flag=bool(self.exhaustive_flag.get()),
                    max_seqs=int(self.max_seqs.get()),
                    evalue=float(self.evalue.get()),
                    alignment_type=int(self.alignment_type.get()),
                    cov=float(self.coverage.get()),
                    cov_mode=int(self.cov_mode.get()),
                    gpu_enabled=bool(self.gpu_enabled.get()),
                    cluster_search=bool(self.cluster_search.get()),
                )
                cmds.append(cmd)
                # Write a per-output-file command log that matches the .m8 filename
                try:
                    per_log = out_file.parent / (out_file.stem + ".cmdlog.txt")
                    write_named_command_log(per_log, cmd)
                except Exception as e:
                    self._log(f"[WARN] Could not write per-file cmdlog: {e}\n")


            self._run_background(cmds, note=f"Running {len(files)} search(es)…", log_dir=None)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _run_background(self, cmd_list: list[list[str]], note: str = "", log_dir: Path | None = None):
        if self.worker and self.worker.is_alive():
            messagebox.showinfo("Busy", "A job is already running.")
            return
        self._log(f"[INFO] {note}\n")
        self.worker = threading.Thread(target=self._worker_fn, args=(cmd_list, log_dir), daemon=True)
        self.worker.start()

    def _worker_fn(self, cmd_list: list[list[str]], log_dir: Path | None):
        try:
            for cmd in cmd_list:
                if log_dir is not None:
                    try:
                        write_command_log(log_dir, cmd)
                    except Exception as e:
                        self.log_queue.put(f"[WARN] Could not write command log: {e}\n")
                run_cmd_streamed(cmd, self.log_queue)
        except Exception as e:
            self.log_queue.put(f"[ERROR] {e}\n")

    def _poll_log_queue(self):
        try:
            while True:
                line = self.log_queue.get_nowait()
                self._log(line)
        except queue.Empty:
            pass
        self.after(100, self._poll_log_queue)

    def _log(self, text: str):
        self.txt.insert("end", text)
        self.txt.see("end")

# -------------------- Main --------------------

def main():
    # ensure foldseek exists (helpful error early)
    try:
        which_or_die("foldseek")
    except FileNotFoundError as e:
        messagebox.showerror("Foldseek not found", str(e))
        return

    app = FoldseekGUI()
    app.mainloop()

if __name__ == "__main__":
    main()