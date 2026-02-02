from __future__ import annotations

import subprocess
import sys
import textwrap
import webbrowser
import csv
import random
import json
import urllib.request
import urllib.error
from pathlib import Path
import importlib.util
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk

try:
    # Optional: drag-and-drop support for files
    from tkinterdnd2 import DND_FILES, TkinterDnD  # type: ignore
except Exception:  # pragma: no cover
    DND_FILES = None  # type: ignore
    TkinterDnD = None  # type: ignore

from tkinter.scrolledtext import ScrolledText

from dna_app import (
    analyze_protein_impact,
    detect_mutations_from_alignment,
    find_orfs,
    mutation_rate,
    needleman_wunsch_global,
    parse_fasta_text,
    read_fasta_file,
    reverse_complement,
    split_codons,
    stats,
    summarize_mutations,
    translate_dna,
    validate_dna,
)
from dna_app.mutations import classify_substitution, Mutation


def _format_wrapped(seq: str, width: int = 80) -> str:
    return "\n".join(textwrap.wrap(seq, width=width)) if seq else ""


def _clean_protein_text(text: str) -> str:
    """Accept FASTA-like text or raw AA sequence, return clean uppercase AA string."""
    if not text:
        return ""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return ""

    # Drop FASTA headers
    if any(ln.startswith(">") for ln in lines):
        lines = [ln for ln in lines if not ln.startswith(">")]

    seq = "".join(lines)
    # Remove whitespace and common formatting characters
    seq = "".join(ch for ch in seq if ch.isalpha() or ch in {"*"})
    return seq.upper()


def _looks_like_protein(seq: str) -> bool:
    if not seq:
        return False
    # Allow common ambiguity letters too.
    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")
    return all(ch in allowed for ch in seq) and len(seq) >= 20


# Built-in visualization examples (PDB IDs) for quick demo.
_VIS_PDB_EXAMPLES: list[tuple[str, str]] = [
    ("Crambin (PDB 1CRN)", "1CRN"),
    ("Ubiquitin (PDB 1UBQ)", "1UBQ"),
    ("Myoglobin (PDB 1MBN)", "1MBN"),
    ("Lysozyme (PDB 1LYZ)", "1LYZ"),
    ("Hemoglobin (PDB 4HHB)", "4HHB"),
]


# If tkinterdnd2 is installed, use TkinterDnD.Tk as the base class.
_BaseTk = TkinterDnD.Tk if TkinterDnD is not None else tk.Tk


class DNAApp(_BaseTk):
    def __init__(self) -> None:
        super().__init__()
        self.title("DNA Sequence Analysis & Mutation Detection")
        self.geometry("1100x750")

        self._last_compare_alignment = None
        self._last_compare_mutations: list[Mutation] = []
        self._last_compare_ref: str = ""
        self._last_compare_alt: str = ""

        self.last_protein_single: str = ""
        self.last_protein_ref: str = ""
        self.last_protein_alt: str = ""

        self._configure_style()
        self._build_ui()

        # Enable drag & drop if tkinterdnd2 is installed.
        self._setup_drag_and_drop()

    def _configure_style(self) -> None:
        try:
            style = ttk.Style(self)
            if "clam" in style.theme_names():
                style.theme_use("clam")

            default_font = ("Segoe UI", 10)
            style.configure("TLabel", font=default_font)
            style.configure("TButton", font=default_font, padding=(10, 6))
            style.configure("TCheckbutton", font=default_font)
            style.configure("TNotebook.Tab", font=("Segoe UI", 10, "bold"), padding=(12, 8))
        except Exception:
            # Keep a safe fallback if a theme/font isn't available.
            pass

    def _open_url(self, url: str) -> None:
        try:
            webbrowser.open(url)
        except Exception as e:
            messagebox.showerror("Open URL failed", f"Could not open browser.\n\nError: {e}")

    def _best_last_protein(self) -> str:
        return self.last_protein_alt or self.last_protein_ref or self.last_protein_single

    def _rcsb_sequence_search(self, protein: str, *, identity_cutoff: float, evalue_cutoff: float, rows: int = 25, include_computed: bool = False) -> list[tuple[str, float]]:
        """Return a list of (PDB_ID, score) from RCSB sequence search."""
        payload: dict = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "evalue_cutoff": float(evalue_cutoff),
                    "identity_cutoff": float(identity_cutoff),
                    "sequence_type": "protein",
                    "value": protein,
                },
            },
            "request_options": {
                "scoring_strategy": "sequence",
                "paginate": {"start": 0, "rows": int(rows)},
                "results_verbosity": "minimal",
            },
            "return_type": "entry",
        }
        if include_computed:
            payload["request_options"]["results_content_type"] = ["experimental", "computational"]

        data = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(
            "https://search.rcsb.org/rcsbsearch/v2/query",
            data=data,
            method="POST",
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )

        try:
            with urllib.request.urlopen(req, timeout=25) as resp:
                body = resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            detail = ""
            try:
                detail = e.read().decode("utf-8", errors="replace")
            except Exception:
                pass
            raise RuntimeError(f"RCSB API error (HTTP {e.code}).\n{detail[:5000]}")
        except Exception as e:
            raise RuntimeError(f"Failed to contact RCSB Search API.\n{e}")

        try:
            obj = json.loads(body)
        except Exception:
            raise RuntimeError("RCSB API returned invalid JSON.")

        hits: list[tuple[str, float]] = []
        for item in obj.get("result_set", []) or []:
            pdb_id = (item.get("identifier") or "").strip().upper()
            score = float(item.get("score") or 0.0)
            if pdb_id:
                hits.append((pdb_id, score))
        return hits

    def _load_example_into(self, widget: ScrolledText, example_path: Path) -> None:
        if not example_path.exists():
            messagebox.showerror("Example not found", f"Missing example file:\n{example_path}")
            return
        try:
            rec = read_fasta_file(str(example_path))
        except Exception as e:
            messagebox.showerror("Failed to read example", f"Could not read:\n{example_path}\n\nError: {e}")
            return

        widget.delete("1.0", "end")
        if rec.header:
            widget.insert("1.0", f">{rec.header}\n")
        widget.insert("end", _format_wrapped(rec.sequence, 80))

    def _setup_drag_and_drop(self) -> None:
        """Register drop targets for FASTA file drag-and-drop (optional)."""
        if DND_FILES is None:
            return

        def register(widget: tk.Widget, handler) -> None:
            try:
                widget.drop_target_register(DND_FILES)  # type: ignore[attr-defined]
                widget.dnd_bind("<<Drop>>", handler)  # type: ignore[attr-defined]
            except Exception:
                # Don't break the app if DnD isn't supported in this environment.
                pass

        # Single analysis input
        if hasattr(self, "single_input"):
            register(self.single_input, self._on_drop_single)

        # Compare inputs
        if hasattr(self, "compare_ref"):
            register(self.compare_ref, self._on_drop_compare_ref)
        if hasattr(self, "compare_alt"):
            register(self.compare_alt, self._on_drop_compare_alt)

        # Visualization protein box (optional)
        if hasattr(self, "vis_protein_text"):
            register(self.vis_protein_text, self._on_drop_protein)

    def _parse_dnd_files(self, data: str) -> list[str]:
        """Parse tkinterdnd2 DND_FILES payload into a list of file paths."""
        if not data:
            return []

        # tkinterdnd2 uses Tcl list semantics; splitlist handles braces/spaces safely.
        try:
            parts = list(self.tk.splitlist(data))
        except Exception:
            parts = [data]

        out: list[str] = []
        for p in parts:
            p = str(p).strip()
            if not p:
                continue
            # Remove file:// prefix if present
            if p.lower().startswith("file://"):
                p = p[7:]
                if p.startswith("/") and len(p) > 3 and p[2] == ":":
                    p = p[1:]
            out.append(p)
        return out

    def _load_fasta_file_into_widget(self, widget: ScrolledText, path: str) -> None:
        try:
            record = read_fasta_file(path)
        except Exception as e:
            messagebox.showerror("Failed to read file", f"Could not read FASTA file:\n{path}\n\nError: {e}")
            return

        widget.delete("1.0", "end")
        if record.header:
            widget.insert("1.0", f">{record.header}\n")
        widget.insert("end", _format_wrapped(record.sequence, 80))

    def _on_drop_single(self, event) -> str:
        files = self._parse_dnd_files(getattr(event, "data", ""))
        if not files:
            return "break"
        path = files[0]
        try:
            record = read_fasta_file(path)
        except Exception as e:
            messagebox.showerror("Failed to read file", f"Could not read FASTA file:\n{path}\n\nError: {e}")
            return "break"
        self.single_header.config(text=f"Header: {record.header}" if record.header else "")
        self.single_input.delete("1.0", "end")
        if record.header:
            self.single_input.insert("1.0", f">{record.header}\n")
        self.single_input.insert("end", _format_wrapped(record.sequence, 80))
        return "break"

    def _on_drop_compare_ref(self, event) -> str:
        files = self._parse_dnd_files(getattr(event, "data", ""))
        if not files:
            return "break"
        self._load_fasta_file_into_widget(self.compare_ref, files[0])
        return "break"

    def _on_drop_compare_alt(self, event) -> str:
        files = self._parse_dnd_files(getattr(event, "data", ""))
        if not files:
            return "break"
        self._load_fasta_file_into_widget(self.compare_alt, files[0])
        return "break"

    def _on_drop_protein(self, event) -> str:
        files = self._parse_dnd_files(getattr(event, "data", ""))
        if not files:
            return "break"
        path = files[0]
        try:
            record = read_fasta_file(path)
        except Exception as e:
            messagebox.showerror("Failed to read file", f"Could not read FASTA file:\n{path}\n\nError: {e}")
            return "break"
        # For protein visualization, we just insert the sequence as-is.
        self.vis_protein_text.delete("1.0", "end")
        self.vis_protein_text.insert("1.0", _format_wrapped(record.sequence, 80))
        return "break"

    def _build_ui(self) -> None:
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.nb = ttk.Notebook(self)
        self.nb.grid(row=0, column=0, sticky="nsew")

        self.tab_single = ttk.Frame(self.nb)
        self.tab_compare = ttk.Frame(self.nb)
        self.tab_visual = ttk.Frame(self.nb)
        self.tab_guide = ttk.Frame(self.nb)

        self.nb.add(self.tab_single, text="Single Analysis")
        self.nb.add(self.tab_compare, text="Compare")
        self.nb.add(self.tab_visual, text="Visualization")
        self.nb.add(self.tab_guide, text="Guide")

        self._build_single_tab(self.tab_single)
        self._build_compare_tab(self.tab_compare)
        self._build_visual_tab(self.tab_visual)
        self._build_guide_tab(self.tab_guide)

    # ----------------------
    # Single analysis
    # ----------------------
    def _build_single_tab(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=1)
        parent.columnconfigure(1, weight=1)
        parent.rowconfigure(1, weight=1)
        parent.rowconfigure(3, weight=1)

        top = ttk.Frame(parent)
        top.grid(row=0, column=0, columnspan=2, sticky="ew", padx=10, pady=10)
        top.columnconfigure(1, weight=1)

        ttk.Label(top, text="Input (FASTA text or raw DNA):").grid(row=0, column=0, sticky="w")
        self.single_header = ttk.Label(top, text="")
        self.single_header.grid(row=0, column=1, sticky="w", padx=(10, 0))

        ttk.Label(top, text="Tip: drag & drop a FASTA file into the box.").grid(row=1, column=0, columnspan=2, sticky="w", pady=(6, 0))

        btns = ttk.Frame(top)
        btns.grid(row=0, column=2, sticky="e")

        ttk.Button(btns, text="Load FASTA file", command=self._single_load_file).grid(row=0, column=0, padx=5)
        ttk.Button(btns, text="Load example", command=self._single_load_example).grid(row=0, column=1, padx=5)
        ttk.Button(btns, text="Analyze", command=self._single_analyze).grid(row=0, column=1, padx=5)
        ttk.Button(btns, text="Clear", command=self._single_clear).grid(row=0, column=2, padx=5)

        self.single_input = ScrolledText(parent, height=12)
        self.single_input.grid(row=1, column=0, columnspan=2, sticky="nsew", padx=10)

        opts = ttk.Frame(parent)
        opts.grid(row=2, column=0, columnspan=2, sticky="ew", padx=10, pady=6)
        opts.columnconfigure(0, weight=1)

        self.single_include_reverse = tk.BooleanVar(value=False)
        ttk.Checkbutton(opts, text="Also scan reverse-complement (6-frame ORFs)", variable=self.single_include_reverse).grid(
            row=0, column=0, sticky="w"
        )

        self.single_output = ScrolledText(parent, height=16, state="disabled")
        self.single_output.grid(row=3, column=0, columnspan=2, sticky="nsew", padx=10, pady=(0, 10))

    def _single_clear(self) -> None:
        self.single_input.delete("1.0", "end")
        self.single_header.config(text="")
        self._set_text(self.single_output, "")

    def _single_load_file(self) -> None:
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA", "*.fasta *.fa *.fas *.txt"), ("All files", "*.*")],
        )
        if not path:
            return
        record = read_fasta_file(path)
        self.single_header.config(text=f"Header: {record.header}" if record.header else "")
        self.single_input.delete("1.0", "end")
        self.single_input.insert("1.0", f">{record.header}\n" if record.header else "")
        self.single_input.insert("end", _format_wrapped(record.sequence, 80))

    def _single_load_example(self) -> None:
        example = Path("examples/single_example.fasta")
        self._load_example_into(self.single_input, example)
        try:
            rec = read_fasta_file(str(example))
            self.single_header.config(text=f"Header: {rec.header}" if rec.header else "")
        except Exception:
            self.single_header.config(text="")

    def _single_analyze(self) -> None:
        text = self.single_input.get("1.0", "end").strip()
        record = parse_fasta_text(text)
        seq = record.sequence
        self.single_header.config(text=f"Header: {record.header}" if record.header else "")

        ok, msg = validate_dna(seq)
        if not ok:
            messagebox.showerror("Invalid sequence", msg)
            return

        s = stats(seq)
        rc = reverse_complement(seq)
        orfs = find_orfs(seq, include_reverse_complement=self.single_include_reverse.get())
        # store a protein for visualization/export (longest ORF if possible)
        self.last_protein_single = max((o.protein for o in orfs), key=len, default="")

        codons_frame0 = split_codons(seq, frame=0)
        start_count = sum(1 for c in codons_frame0 if c == "ATG")
        stop_count = sum(1 for c in codons_frame0 if c in {"TAA", "TAG", "TGA"})

        # Translate full sequence in the 3 frames (educational output)
        prot0 = translate_dna(seq, frame=0, stop_at_stop=False)
        prot1 = translate_dna(seq, frame=1, stop_at_stop=False)
        prot2 = translate_dna(seq, frame=2, stop_at_stop=False)

        lines: list[str] = []
        lines.append("PART 1 — Single DNA Sequence Analysis")
        lines.append("=" * 60)
        lines.append(f"Length: {s.length} nt")
        lines.append(f"A:{s.count_a}  T:{s.count_t}  C:{s.count_c}  G:{s.count_g}")
        lines.append(f"GC%: {s.gc_percent:.2f}   AT%: {s.at_percent:.2f}")
        lines.append("")

        lines.append("Reverse-complement:")
        lines.append(_format_wrapped(rc, 80))
        lines.append("")

        lines.append("Codon analysis (frame 1 / offset 0):")
        lines.append(f"Codons: {len(codons_frame0)}")
        lines.append(f"Start codon (ATG) count: {start_count}")
        lines.append(f"Stop codons (TAA/TAG/TGA) count: {stop_count}")
        lines.append("")

        lines.append(f"ORFs detected: {len(orfs)}")
        for idx, o in enumerate(orfs[:10], 1):
            lines.append(
                f"  ORF #{idx}: strand {o.strand}, frame {o.frame+1}, {o.start+1}..{o.end} ({len(o.dna)} nt, {len(o.protein)} aa)"
            )
        if len(orfs) > 10:
            lines.append(f"  ... ({len(orfs)-10} more)")
        lines.append("")

        lines.append("Protein translation (3 frames):")
        lines.append("Frame 1:")
        lines.append(_format_wrapped(prot0, 80))
        lines.append("Frame 2:")
        lines.append(_format_wrapped(prot1, 80))
        lines.append("Frame 3:")
        lines.append(_format_wrapped(prot2, 80))

        self._set_text(self.single_output, "\n".join(lines))

    # ----------------------
    # Compare
    # ----------------------
    def _build_compare_tab(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=1)
        parent.columnconfigure(1, weight=1)
        parent.rowconfigure(1, weight=1)
        parent.rowconfigure(3, weight=2)

        top = ttk.Frame(parent)
        top.grid(row=0, column=0, columnspan=2, sticky="ew", padx=10, pady=10)
        top.columnconfigure(0, weight=1)

        ttk.Label(top, text="Tip: drag & drop FASTA files into REF/MUT boxes.").grid(row=1, column=0, sticky="w", pady=(6, 0))

        btns = ttk.Frame(top)
        btns.grid(row=0, column=1, sticky="e")
        ttk.Button(btns, text="Load REF FASTA", command=lambda: self._compare_load(which="ref")).grid(row=0, column=0, padx=5)
        ttk.Button(btns, text="Load MUT FASTA", command=lambda: self._compare_load(which="alt")).grid(row=0, column=1, padx=5)
        ttk.Button(btns, text="Load random example", command=self._compare_load_examples).grid(row=0, column=2, padx=5)
        ttk.Button(btns, text="Compare", command=self._compare_run).grid(row=0, column=3, padx=5)
        ttk.Button(btns, text="Mutation table", command=self._show_mutation_table).grid(row=0, column=4, padx=5)
        ttk.Button(btns, text="Export CSV", command=self._export_mutations_csv).grid(row=0, column=5, padx=5)
        ttk.Button(btns, text="Copy report", command=self._copy_compare_report).grid(row=0, column=6, padx=5)
        ttk.Button(btns, text="Clear", command=self._compare_clear).grid(row=0, column=7, padx=5)

        ttk.Label(parent, text="Reference sequence (FASTA or raw):").grid(row=1, column=0, sticky="w", padx=10)
        ttk.Label(parent, text="Mutated sequence (FASTA or raw):").grid(row=1, column=1, sticky="w", padx=10)

        self.compare_ref = ScrolledText(parent, height=10)
        self.compare_alt = ScrolledText(parent, height=10)
        self.compare_ref.grid(row=2, column=0, sticky="nsew", padx=10)
        self.compare_alt.grid(row=2, column=1, sticky="nsew", padx=10)

        legend = ttk.Frame(parent)
        legend.grid(row=3, column=0, columnspan=2, sticky="ew", padx=10, pady=(6, 4))
        legend.columnconfigure(10, weight=1)

        ttk.Label(legend, text="Legend:").grid(row=0, column=0, sticky="w", padx=(0, 8))
        self._legend_chip(legend, "Match", "#166534", 1)
        self._legend_chip(legend, "Transition", "#1d4ed8", 2)
        self._legend_chip(legend, "Transversion", "#6d28d9", 3)
        self._legend_chip(legend, "Insertion", "#b45309", 4)
        self._legend_chip(legend, "Deletion", "#b91c1c", 5)
        self._legend_chip(legend, "Gap", "#6b7280", 6)
        ttk.Label(legend, text="Tip: scroll alignment; use Mutation table for exact positions.").grid(
            row=0, column=10, sticky="e"
        )

        self.compare_output = ScrolledText(parent, height=20, state="disabled")
        self.compare_output.grid(row=4, column=0, columnspan=2, sticky="nsew", padx=10, pady=(0, 10))

    def _legend_chip(self, parent: ttk.Frame, text: str, color: str, col: int) -> None:
        chip = tk.Label(parent, text=f"  {text}  ", bg=color, fg="white")
        chip.grid(row=0, column=col, sticky="w", padx=4)

    def _compare_clear(self) -> None:
        self.compare_ref.delete("1.0", "end")
        self.compare_alt.delete("1.0", "end")
        self._set_text(self.compare_output, "")
        self._last_compare_alignment = None
        self._last_compare_mutations = []
        self._last_compare_ref = ""
        self._last_compare_alt = ""

    def _compare_load_examples(self) -> None:
        # Load a random example pair from known datasets.
        candidates: list[tuple[Path, Path]] = []
        known_pairs = [
            (Path("examples/reference.fasta"), Path("examples/mutated.fasta")),
            (Path("TP4/origin.fasta"), Path("TP4/mutation.fasta")),
        ]
        for rp, ap in known_pairs:
            if rp.exists() and ap.exists():
                candidates.append((rp, ap))

        # Also auto-discover pairs in examples/ by naming convention: reference*.fasta <-> mutated*.fasta
        ex_dir = Path("examples")
        if ex_dir.exists():
            refs = list(ex_dir.glob("reference*.fasta")) + list(ex_dir.glob("reference*.fa"))
            muts = list(ex_dir.glob("mutated*.fasta")) + list(ex_dir.glob("mutated*.fa"))
            ref_map: dict[str, Path] = {p.stem.replace("reference", "", 1): p for p in refs}
            mut_map: dict[str, Path] = {p.stem.replace("mutated", "", 1): p for p in muts}
            for key, rp in ref_map.items():
                ap = mut_map.get(key)
                if ap is not None:
                    candidates.append((rp, ap))

        if not candidates:
            messagebox.showerror(
                "No examples found",
                "No example pairs were found.\n\nExpected one of:\n"
                "- examples/reference.fasta + examples/mutated.fasta\n"
                "- TP4/origin.fasta + TP4/mutation.fasta",
            )
            return

        ref_path, alt_path = random.choice(candidates)
        self._load_example_into(self.compare_ref, ref_path)
        self._load_example_into(self.compare_alt, alt_path)
        messagebox.showinfo("Example loaded", f"Loaded:\nREF: {ref_path}\nMUT: {alt_path}")

    def _compare_load(self, which: str) -> None:
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA", "*.fasta *.fa *.fas *.txt"), ("All files", "*.*")],
        )
        if not path:
            return
        record = read_fasta_file(path)
        box = self.compare_ref if which == "ref" else self.compare_alt
        box.delete("1.0", "end")
        if record.header:
            box.insert("1.0", f">{record.header}\n")
        box.insert("end", _format_wrapped(record.sequence, 80))

    def _compare_run(self) -> None:
        ref = parse_fasta_text(self.compare_ref.get("1.0", "end")).sequence
        alt = parse_fasta_text(self.compare_alt.get("1.0", "end")).sequence

        ok1, msg1 = validate_dna(ref)
        ok2, msg2 = validate_dna(alt)
        if not ok1:
            messagebox.showerror("Invalid reference", msg1)
            return
        if not ok2:
            messagebox.showerror("Invalid mutated sequence", msg2)
            return

        aln = needleman_wunsch_global(ref, alt)
        muts = detect_mutations_from_alignment(aln)
        summary = summarize_mutations(muts)
        rate = mutation_rate(muts, ref_length=len(ref))

        impact = analyze_protein_impact(ref, alt)

        # store proteins for visualization/export
        self.last_protein_ref = impact.ref_protein
        self.last_protein_alt = impact.alt_protein

        # store compare artifacts for mutation table/export
        self._last_compare_alignment = aln
        self._last_compare_mutations = muts
        self._last_compare_ref = ref
        self._last_compare_alt = alt

        # ORF comparison (counts)
        orfs_ref = find_orfs(ref, include_reverse_complement=False)
        orfs_alt = find_orfs(alt, include_reverse_complement=False)

        self.compare_output.config(state="normal")
        self.compare_output.delete("1.0", "end")

        self.compare_output.tag_configure("mono", font=("Consolas", 10))

        # per-base coloring tags
        self.compare_output.tag_configure("match", foreground="#166534")
        self.compare_output.tag_configure("transition", foreground="#1d4ed8")
        self.compare_output.tag_configure("transversion", foreground="#6d28d9")
        self.compare_output.tag_configure("insertion", foreground="#b45309")
        self.compare_output.tag_configure("deletion", foreground="#b91c1c")
        self.compare_output.tag_configure("gap", foreground="#6b7280")
        self.compare_output.tag_configure("mut_col", background="#fff1f2")

        lines: list[str] = []
        lines.append("PART 2 — Two DNA Sequence Comparison")
        lines.append("=" * 60)
        lines.append(f"Reference length: {len(ref)} nt")
        lines.append(f"Mutated length:   {len(alt)} nt")
        lines.append("")

        lines.append("Mutations summary:")
        lines.append(f"  Total: {len(muts)}")
        lines.append(
            f"  Substitutions: {summary['substitution']} (transition {summary['transition']}, transversion {summary['transversion']})"
        )
        lines.append(f"  Insertions: {summary['insertion']}")
        lines.append(f"  Deletions: {summary['deletion']}")
        lines.append(f"  Mutation rate (per ref base): {rate:.4f}")
        lines.append("")

        lines.append("ORF comparison (+ strand):")
        lines.append(f"  ORFs (ref): {len(orfs_ref)}")
        lines.append(f"  ORFs (mut): {len(orfs_alt)}")
        lines.append("")

        lines.append("Protein impact (heuristic, longest ORF):")
        lines.append(f"  Impact: {impact.impact}")
        lines.append(f"  Notes: {impact.notes}")
        if impact.ref_orf_len_nt is not None and impact.alt_orf_len_nt is not None:
            lines.append(f"  ORF length (ref → mut): {impact.ref_orf_len_nt} nt → {impact.alt_orf_len_nt} nt")
        lines.append("")

        lines.append("PART 3 — Alignment visualization")
        lines.append("=" * 60)
        lines.append("Legend: matches in green; transitions blue; transversions purple; insertions orange; deletions red.\n")

        self.compare_output.insert("end", "\n".join(lines) + "\n")

        # Precompute per-column reference/alt positions and mutation types for rich coloring.
        ref_pos_at: list[int | None] = []
        alt_pos_at: list[int | None] = []
        kind_at: list[str] = []

        ref_i = 0
        alt_i = 0
        for rch, ach in zip(aln.aligned_ref, aln.aligned_alt, strict=True):
            if rch != "-":
                ref_i += 1
                rp: int | None = ref_i
            else:
                rp = None

            if ach != "-":
                alt_i += 1
                ap: int | None = alt_i
            else:
                ap = None

            ref_pos_at.append(rp)
            alt_pos_at.append(ap)

            if rch == ach:
                kind_at.append("match")
            elif rch == "-":
                kind_at.append("insertion")
            elif ach == "-":
                kind_at.append("deletion")
            else:
                cls = classify_substitution(rch, ach)
                kind_at.append(cls if cls in {"transition", "transversion"} else "transversion")

        # chunk alignment in blocks for readability, and color bases directly
        block = 80
        for k in range(0, len(aln.aligned_ref), block):
            r = aln.aligned_ref[k : k + block]
            m = aln.midline[k : k + block]
            a = aln.aligned_alt[k : k + block]

            rp = ref_pos_at[k : k + block]
            ap = alt_pos_at[k : k + block]
            kp = kind_at[k : k + block]

            ref_first = next((x for x in rp if x is not None), None)
            ref_last = next((x for x in reversed(rp) if x is not None), None)
            alt_first = next((x for x in ap if x is not None), None)
            alt_last = next((x for x in reversed(ap) if x is not None), None)

            hdr = f"Cols {k+1}-{k+len(r)}"
            if ref_first is not None and ref_last is not None:
                hdr += f"  |  REF {ref_first}-{ref_last}"
            if alt_first is not None and alt_last is not None:
                hdr += f"  |  MUT {alt_first}-{alt_last}"

            self.compare_output.insert("end", hdr + "\n", ("mono",))

            self.compare_output.insert("end", "REF ", ("mono",))
            self._insert_colored_bases(self.compare_output, r, kp, for_row="ref")
            self.compare_output.insert("end", "\n", ("mono",))

            self.compare_output.insert("end", "MID ", ("mono",))
            for idx, ch in enumerate(m):
                tag = "match" if ch == "|" else ("deletion" if ch == "*" else "gap")
                extra = ("mut_col",) if ch == "*" else ()
                self.compare_output.insert("end", ch, ("mono", tag, *extra))
            self.compare_output.insert("end", "\n", ("mono",))

            self.compare_output.insert("end", "MUT ", ("mono",))
            self._insert_colored_bases(self.compare_output, a, kp, for_row="alt")
            self.compare_output.insert("end", "\n\n", ("mono",))

        if muts:
            self.compare_output.insert("end", "Mutations (first 50):\n")
            for mu in muts[:50]:
                if mu.kind == "substitution":
                    self.compare_output.insert(
                        "end",
                        f"  Substitution ref:{mu.ref_pos} {mu.ref_base}→{mu.alt_base} ({mu.substitution_class})\n",
                    )
                elif mu.kind == "insertion":
                    self.compare_output.insert("end", f"  Insertion after ref:{mu.ref_pos} +{mu.alt_base} (alt pos {mu.alt_pos})\n")
                else:
                    self.compare_output.insert("end", f"  Deletion ref:{mu.ref_pos} -{mu.ref_base} (alt pos {mu.alt_pos})\n")
            if len(muts) > 50:
                self.compare_output.insert("end", f"  ... ({len(muts)-50} more)\n")

        self.compare_output.insert("end", "\nProtein comparison (longest ORF):\nREF:\n")
        self.compare_output.insert("end", _format_wrapped(impact.ref_protein, 80) + "\n")
        self.compare_output.insert("end", "MUT:\n")
        self.compare_output.insert("end", _format_wrapped(impact.alt_protein, 80) + "\n")

        self.compare_output.config(state="disabled")

    def _insert_colored_bases(self, widget: ScrolledText, bases: str, kind_at: list[str], *, for_row: str) -> None:
        for ch, kind in zip(bases, kind_at, strict=True):
            tags = ["mono"]
            if ch == "-":
                tags.append("gap")
            elif kind == "match":
                tags.append("match")
            elif kind == "insertion":
                tags.append("insertion")
            elif kind == "deletion":
                tags.append("deletion")
            elif kind == "transition":
                tags.append("transition")
            else:
                tags.append("transversion")

            # Light background highlight on mutated columns (except pure gaps)
            if kind != "match" and ch != "-":
                tags.append("mut_col")

            widget.insert("end", ch, tuple(tags))

    def _show_mutation_table(self) -> None:
        if not self._last_compare_mutations:
            messagebox.showinfo("No mutations", "Run Compare first to generate mutations.")
            return

        win = tk.Toplevel(self)
        win.title("Mutations")
        win.geometry("900x500")

        cols = ("kind", "ref_pos", "alt_pos", "ref", "alt", "class")
        tree = ttk.Treeview(win, columns=cols, show="headings")
        tree.heading("kind", text="Type")
        tree.heading("ref_pos", text="Ref pos")
        tree.heading("alt_pos", text="Mut pos")
        tree.heading("ref", text="Ref")
        tree.heading("alt", text="Mut")
        tree.heading("class", text="Class")

        tree.column("kind", width=120, anchor="w")
        tree.column("ref_pos", width=90, anchor="e")
        tree.column("alt_pos", width=90, anchor="e")
        tree.column("ref", width=70, anchor="center")
        tree.column("alt", width=70, anchor="center")
        tree.column("class", width=120, anchor="w")

        yscroll = ttk.Scrollbar(win, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=yscroll.set)

        tree.grid(row=0, column=0, sticky="nsew")
        yscroll.grid(row=0, column=1, sticky="ns")

        win.columnconfigure(0, weight=1)
        win.rowconfigure(0, weight=1)

        for mu in self._last_compare_mutations:
            cls = mu.substitution_class if mu.kind == "substitution" else ""
            tree.insert(
                "",
                "end",
                values=(mu.kind, mu.ref_pos if mu.ref_pos is not None else "", mu.alt_pos if mu.alt_pos is not None else "", mu.ref_base, mu.alt_base, cls),
            )

        controls = ttk.Frame(win)
        controls.grid(row=1, column=0, columnspan=2, sticky="ew", padx=10, pady=10)
        controls.columnconfigure(0, weight=1)

        ttk.Label(controls, text=f"Total mutations: {len(self._last_compare_mutations)}").grid(row=0, column=0, sticky="w")
        ttk.Button(controls, text="Export CSV", command=self._export_mutations_csv).grid(row=0, column=1, padx=5)
        ttk.Button(controls, text="Copy as TSV", command=lambda: self._copy_mutations_tsv()).grid(row=0, column=2, padx=5)

    def _copy_mutations_tsv(self) -> None:
        if not self._last_compare_mutations:
            return
        rows = ["kind\tref_pos\talt_pos\tref\talt\tclass"]
        for mu in self._last_compare_mutations:
            cls = mu.substitution_class if mu.kind == "substitution" else ""
            rows.append(
                "\t".join(
                    [
                        mu.kind,
                        str(mu.ref_pos or ""),
                        str(mu.alt_pos or ""),
                        mu.ref_base,
                        mu.alt_base,
                        cls,
                    ]
                )
            )
        text = "\n".join(rows)
        self.clipboard_clear()
        self.clipboard_append(text)
        self.update()

    def _export_mutations_csv(self) -> None:
        if not self._last_compare_mutations:
            messagebox.showinfo("No mutations", "Run Compare first to generate mutations.")
            return
        path = filedialog.asksaveasfilename(
            title="Export mutations as CSV",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["kind", "ref_pos", "alt_pos", "ref_base", "alt_base", "substitution_class"])
                for mu in self._last_compare_mutations:
                    w.writerow([
                        mu.kind,
                        mu.ref_pos if mu.ref_pos is not None else "",
                        mu.alt_pos if mu.alt_pos is not None else "",
                        mu.ref_base,
                        mu.alt_base,
                        mu.substitution_class if mu.kind == "substitution" else "",
                    ])
        except Exception as e:
            messagebox.showerror("Export failed", f"Could not export CSV.\n\nError: {e}")
            return

        messagebox.showinfo("Exported", f"Saved: {path}")

    def _copy_compare_report(self) -> None:
        text = self.compare_output.get("1.0", "end").strip()
        if not text:
            return
        self.clipboard_clear()
        self.clipboard_append(text)
        self.update()

    # ----------------------
    # Visualization tab (advanced hook)
    # ----------------------
    def _build_visual_tab(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(1, weight=0)
        parent.rowconfigure(2, weight=1)

        hdr = ttk.Frame(parent)
        hdr.grid(row=0, column=0, sticky="ew", padx=10, pady=10)
        hdr.columnconfigure(0, weight=1)

        ttk.Label(
            hdr,
            text=(
                "3D Protein Visualization (Advanced):\n"
                "This project can integrate external viewers (Mol*, NGL, RCSB, AlphaFold DB).\n"
                "In this lightweight version, we open web viewers where you can paste a protein/PDB ID.\n"
                "Tip: open the Guide tab for step-by-step instructions."
            ),
            justify="left",
        ).grid(row=0, column=0, sticky="w")

        btns = ttk.Frame(hdr)
        btns.grid(row=0, column=1, sticky="e")
        ttk.Button(btns, text="Open Mol* Viewer", command=lambda: webbrowser.open("https://molstar.org/viewer/")).grid(
            row=0, column=0, padx=5
        )
        ttk.Button(btns, text="Open RCSB PDB", command=lambda: webbrowser.open("https://www.rcsb.org/")).grid(row=0, column=1, padx=5)
        ttk.Button(btns, text="Open AlphaFold DB", command=lambda: webbrowser.open("https://alphafold.ebi.ac.uk/")).grid(row=0, column=2, padx=5)
        ttk.Button(btns, text="Open Guide", command=self._go_to_guide).grid(row=0, column=3, padx=5)

        tools = ttk.Frame(parent)
        tools.grid(row=1, column=0, sticky="ew", padx=10, pady=(0, 10))
        tools.columnconfigure(1, weight=1)

        ttk.Label(tools, text="PDB ID (optional):").grid(row=0, column=0, sticky="w")
        self.pdb_id = ttk.Entry(tools)
        self.pdb_id.grid(row=0, column=1, sticky="ew", padx=8)
        ttk.Button(tools, text="Open embedded 3D", command=self._open_embedded_3d).grid(row=0, column=2, padx=5)

        ttk.Button(tools, text="Save last protein FASTA", command=self._save_last_protein_fasta).grid(
            row=0, column=3, padx=5
        )

        ttk.Separator(tools, orient="horizontal").grid(row=1, column=0, columnspan=4, sticky="ew", pady=8)

        # Quick examples (PDB IDs)
        examples = ttk.Frame(tools)
        examples.grid(row=2, column=0, columnspan=4, sticky="ew")
        examples.columnconfigure(1, weight=1)

        ttk.Label(examples, text="Examples (PDB IDs):").grid(row=0, column=0, sticky="w")
        self.vis_example = tk.StringVar(value=_VIS_PDB_EXAMPLES[0][0] if _VIS_PDB_EXAMPLES else "")
        self.vis_example_combo = ttk.Combobox(
            examples,
            textvariable=self.vis_example,
            state="readonly",
            values=[label for label, _pdb in _VIS_PDB_EXAMPLES],
        )
        self.vis_example_combo.grid(row=0, column=1, sticky="ew", padx=8)
        ttk.Button(examples, text="Use example", command=self._vis_use_example).grid(row=0, column=2, padx=5)
        ttk.Button(examples, text="Open in browser", command=self._vis_open_example_browser).grid(row=0, column=3, padx=5)

        # --- Sequence → structure helper (RCSB sequence search) ---
        seqbox = ttk.Labelframe(parent, text="Find a structure automatically (sequence similarity search)")
        seqbox.grid(row=2, column=0, sticky="nsew", padx=10, pady=(0, 10))
        seqbox.columnconfigure(0, weight=1)
        seqbox.rowconfigure(1, weight=1)

        seq_top = ttk.Frame(seqbox)
        seq_top.grid(row=0, column=0, sticky="ew", padx=10, pady=8)
        seq_top.columnconfigure(1, weight=1)

        ttk.Label(seq_top, text="Protein source:").grid(row=0, column=0, sticky="w")
        self.vis_protein_source = tk.StringVar(value="auto")
        src = ttk.Frame(seq_top)
        src.grid(row=0, column=1, sticky="w")
        ttk.Radiobutton(src, text="Auto (best)", variable=self.vis_protein_source, value="auto").grid(row=0, column=0, padx=(0, 8))
        ttk.Radiobutton(src, text="Single", variable=self.vis_protein_source, value="single").grid(row=0, column=1, padx=(0, 8))
        ttk.Radiobutton(src, text="Compare REF", variable=self.vis_protein_source, value="ref").grid(row=0, column=2, padx=(0, 8))
        ttk.Radiobutton(src, text="Compare MUT", variable=self.vis_protein_source, value="alt").grid(row=0, column=3)

        seq_btns = ttk.Frame(seq_top)
        seq_btns.grid(row=0, column=2, sticky="e")
        ttk.Button(seq_btns, text="Fill from last", command=self._vis_fill_protein).grid(row=0, column=0, padx=5)
        ttk.Button(seq_btns, text="Copy protein", command=self._vis_copy_protein).grid(row=0, column=1, padx=5)

        self.vis_protein_text = ScrolledText(seqbox, height=6)
        self.vis_protein_text.grid(row=1, column=0, sticky="nsew", padx=10, pady=(0, 8))

        params = ttk.Frame(seqbox)
        params.grid(row=2, column=0, sticky="ew", padx=10, pady=(0, 8))
        params.columnconfigure(9, weight=1)

        ttk.Label(params, text="Identity cutoff:").grid(row=0, column=0, sticky="w")
        self.vis_identity = tk.DoubleVar(value=0.3)
        ttk.Entry(params, textvariable=self.vis_identity, width=6).grid(row=0, column=1, padx=(6, 12), sticky="w")

        ttk.Label(params, text="E-value cutoff:").grid(row=0, column=2, sticky="w")
        self.vis_evalue = tk.DoubleVar(value=0.1)
        ttk.Entry(params, textvariable=self.vis_evalue, width=8).grid(row=0, column=3, padx=(6, 12), sticky="w")

        ttk.Label(params, text="Max hits:").grid(row=0, column=4, sticky="w")
        self.vis_rows = tk.IntVar(value=25)
        ttk.Entry(params, textvariable=self.vis_rows, width=6).grid(row=0, column=5, padx=(6, 12), sticky="w")

        self.vis_include_computed = tk.BooleanVar(value=False)
        ttk.Checkbutton(params, text="Include computed models", variable=self.vis_include_computed).grid(row=0, column=6, sticky="w")

        ttk.Button(params, text="Search RCSB", command=self._vis_search_rcsb).grid(row=0, column=7, padx=6)
        ttk.Button(params, text="Open selected in Mol*", command=self._vis_open_selected_molstar).grid(row=0, column=8, padx=6)
        ttk.Button(params, text="Open selected in RCSB", command=self._vis_open_selected_rcsb).grid(row=0, column=9, sticky="e")

        results = ttk.Frame(seqbox)
        results.grid(row=3, column=0, sticky="nsew", padx=10, pady=(0, 10))
        results.columnconfigure(0, weight=1)
        results.rowconfigure(0, weight=1)

        self.vis_hits_list = tk.Listbox(results, height=7)
        self.vis_hits_list.grid(row=0, column=0, sticky="nsew")
        yscroll = ttk.Scrollbar(results, orient="vertical", command=self.vis_hits_list.yview)
        yscroll.grid(row=0, column=1, sticky="ns")
        self.vis_hits_list.configure(yscrollcommand=yscroll.set)
        self.vis_hits_list.bind("<Double-Button-1>", lambda _e: self._vis_open_selected_molstar())

        hint = ttk.Label(
            seqbox,
            text=(
                "Note: Mol* needs a structure (PDB/mmCIF). This button searches RCSB for similar protein sequences, "
                "then you can open a matching PDB entry directly in Mol*."
            ),
            wraplength=980,
            justify="left",
        )
        hint.grid(row=4, column=0, sticky="w", padx=10, pady=(0, 10))

        info = ttk.Label(
            parent,
            text=(
                "Embedded 3D uses Mol* inside the app via WebView. "
                "If it fails, install pywebview (pip install pywebview) or use the browser buttons above."
            ),
            wraplength=1000,
            justify="left",
        )
        info.grid(row=3, column=0, sticky="w", padx=10, pady=(0, 10))

        # Prefill protein text if available
        self._vis_fill_protein(silent=True)

    def _vis_use_example(self) -> None:
        label = (self.vis_example.get() if hasattr(self, "vis_example") else "").strip()
        pdb = next((p for (lbl, p) in _VIS_PDB_EXAMPLES if lbl == label), "")
        if not pdb:
            return
        try:
            self.pdb_id.delete(0, "end")
            self.pdb_id.insert(0, pdb)
        except Exception:
            pass

        # Also select it in the hits list if present
        if hasattr(self, "vis_hits_list"):
            try:
                self.vis_hits_list.delete(0, "end")
                self.vis_hits_list.insert("end", f"{pdb}\t(example)")
                self.vis_hits_list.selection_clear(0, "end")
                self.vis_hits_list.selection_set(0)
            except Exception:
                pass

    def _vis_open_example_browser(self) -> None:
        label = (self.vis_example.get() if hasattr(self, "vis_example") else "").strip()
        pdb = next((p for (lbl, p) in _VIS_PDB_EXAMPLES if lbl == label), "")
        if not pdb:
            return
        self._open_url(f"https://molstar.org/viewer/?pdb={pdb}")

    def _vis_fill_protein(self, silent: bool = False) -> None:
        src = (self.vis_protein_source.get() if hasattr(self, "vis_protein_source") else "auto").strip()
        if src == "single":
            protein = self.last_protein_single
        elif src == "ref":
            protein = self.last_protein_ref
        elif src == "alt":
            protein = self.last_protein_alt
        else:
            protein = self._best_last_protein()

        if not protein:
            if not silent:
                messagebox.showinfo(
                    "No protein available",
                    "Run Single Analysis or Compare first so the app can extract a protein from an ORF.\n\n"
                    "Or paste a protein sequence manually in the box.",
                )
            return

        self.vis_protein_text.delete("1.0", "end")
        self.vis_protein_text.insert("1.0", _format_wrapped(protein, 80))

    def _vis_copy_protein(self) -> None:
        seq = _clean_protein_text(self.vis_protein_text.get("1.0", "end"))
        if not _looks_like_protein(seq):
            messagebox.showinfo("Nothing to copy", "Protein box is empty (or too short).")
            return
        self.clipboard_clear()
        self.clipboard_append(seq)
        self.update()

    def _vis_search_rcsb(self) -> None:
        protein = _clean_protein_text(self.vis_protein_text.get("1.0", "end"))
        if not _looks_like_protein(protein):
            messagebox.showerror(
                "Invalid protein",
                "Provide a protein amino-acid sequence (FASTA or raw).\n"
                "It must be at least ~20 aa and contain only amino-acid letters.",
            )
            return

        try:
            identity = float(self.vis_identity.get())
            evalue = float(self.vis_evalue.get())
            rows = int(self.vis_rows.get())
        except Exception:
            messagebox.showerror("Invalid parameters", "Identity/E-value/Max hits must be numbers.")
            return

        if not (0.0 < identity <= 1.0):
            messagebox.showerror("Invalid identity cutoff", "Identity cutoff must be between 0 and 1.")
            return
        if evalue <= 0:
            messagebox.showerror("Invalid E-value cutoff", "E-value cutoff must be > 0.")
            return
        if rows <= 0 or rows > 100:
            messagebox.showerror("Invalid max hits", "Max hits must be between 1 and 100.")
            return

        self.vis_hits_list.delete(0, "end")
        self.vis_hits_list.insert("end", "Searching RCSB…")
        self.update_idletasks()

        try:
            hits = self._rcsb_sequence_search(
                protein,
                identity_cutoff=identity,
                evalue_cutoff=evalue,
                rows=rows,
                include_computed=bool(self.vis_include_computed.get()),
            )
        except Exception as e:
            self.vis_hits_list.delete(0, "end")
            messagebox.showerror("RCSB search failed", str(e))
            return

        self.vis_hits_list.delete(0, "end")
        if not hits:
            self.vis_hits_list.insert("end", "No hits. Try lowering identity cutoff or increasing E-value.")
            return
        for pdb_id, score in hits:
            self.vis_hits_list.insert("end", f"{pdb_id}\t(score {score:.3f})")

    def _vis_get_selected_pdb_id(self) -> str:
        try:
            sel = self.vis_hits_list.curselection()
            if not sel:
                return ""
            line = self.vis_hits_list.get(sel[0])
        except Exception:
            return ""

        pdb = (line.split("\t", 1)[0] or "").strip().upper()
        if len(pdb) == 4 and pdb.isalnum():
            return pdb
        return ""

    def _vis_open_selected_molstar(self) -> None:
        pdb = self._vis_get_selected_pdb_id()
        if not pdb:
            messagebox.showinfo("No selection", "Select a hit first (then double-click or click Open).")
            return
        self._open_url(f"https://molstar.org/viewer/?pdb={pdb}")

    def _vis_open_selected_rcsb(self) -> None:
        pdb = self._vis_get_selected_pdb_id()
        if not pdb:
            messagebox.showinfo("No selection", "Select a hit first.")
            return
        self._open_url(f"https://www.rcsb.org/structure/{pdb}")

    # ----------------------
    # Guide tab
    # ----------------------
    def _build_guide_tab(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(1, weight=1)

        top = ttk.Frame(parent)
        top.grid(row=0, column=0, sticky="ew", padx=10, pady=10)
        top.columnconfigure(0, weight=1)

        ttk.Label(
            top,
            text=(
                "Quick Guide — How to use the app + online tools\n"
                "This tab explains how to read the colored alignment and how to use Mol*/RCSB/AlphaFold with your results."
            ),
            justify="left",
        ).grid(row=0, column=0, sticky="w")

        links = ttk.Frame(top)
        links.grid(row=0, column=1, sticky="e")
        ttk.Button(links, text="Mol*", command=lambda: self._open_url("https://molstar.org/viewer/")).grid(row=0, column=0, padx=5)
        ttk.Button(links, text="RCSB PDB", command=lambda: self._open_url("https://www.rcsb.org/")).grid(row=0, column=1, padx=5)
        ttk.Button(links, text="AlphaFold", command=lambda: self._open_url("https://alphafold.ebi.ac.uk/")).grid(row=0, column=2, padx=5)

        guide = ScrolledText(parent, height=22)
        guide.grid(row=1, column=0, sticky="nsew", padx=10, pady=(0, 10))
        guide.configure(font=("Segoe UI", 10))

        content = []
        content.append("1) Compare two sequences")
        content.append("- Go to the Compare tab")
        content.append("- Load REF and MUT FASTA (or click 'Load random example')")
        content.append("- You can also drag & drop FASTA files into the input boxes")
        content.append("- Click Compare")
        content.append("")
        content.append("2) Read the colored alignment")
        content.append("- Green: exact match")
        content.append("- Blue: transition (A↔G or C↔T)")
        content.append("- Purple: transversion (all other substitutions)")
        content.append("- Orange: insertion (extra base in MUT relative to REF)")
        content.append("- Red: deletion (missing base in MUT relative to REF)")
        content.append("- Gray '-': a gap introduced by alignment")
        content.append("")
        content.append("Tips:")
        content.append("- Use 'Mutation table' for exact positions and types")
        content.append("- Use 'Export CSV' to save the mutation list")
        content.append("- Use 'Copy report' to paste into your report")
        content.append("")
        content.append("3) Visualization — Mol* / RCSB / AlphaFold (step-by-step)")
        content.append("A) If you have a PDB ID (best case)")
        content.append("- Go to Visualization tab")
        content.append("- Paste the PDB ID (example: 1CRN) and click 'Open embedded 3D'")
        content.append("- Or click 'Open Mol* Viewer' to open it in your browser")
        content.append("")
        content.append("B) If you do NOT have a PDB ID")
        content.append("- Go to Visualization tab → use 'Find a structure automatically' (sequence similarity search)")
        content.append("  - Click 'Fill from last' to use the protein extracted from the longest ORF")
        content.append("  - Click 'Search RCSB' to get PDB IDs")
        content.append("  - Select a hit and click 'Open selected in Mol*'")
        content.append("- Or use 'Save last protein FASTA' to export the protein and search manually")
        content.append("- AlphaFold DB works when a predicted structure already exists for your protein")
        content.append("")
        content.append("C) How to color/highlight in Mol*")
        content.append("- In Mol*, open the 'Components' panel and select the polymer chain")
        content.append("- Use 'Color Theme' (e.g., by secondary structure, by chain, by residue properties)")
        content.append("- To highlight residues: use selection (by residue range) then apply a representation/color")
        content.append("")
        content.append("D) Practical tip")
        content.append("- This app highlights nucleotide mutations in the alignment; mapping them to 3D requires mapping DNA→protein residue positions")
        content.append("  (future feature: automatically map ORF codon index to residue index for highlighting in Mol*)")
        content.append("")
        content.append("4) Workflow recommendation")
        content.append("- Start with Single Analysis to validate/understand a sequence")
        content.append("- Then Compare for mutation detection and a clean summary")
        content.append("- Finally open Mol*/RCSB/AlphaFold if you want a structure-level view")

        guide.insert("1.0", "\n".join(content))
        guide.configure(state="disabled")

    def _go_to_guide(self) -> None:
        try:
            self.nb.select(self.tab_guide)
        except Exception:
            pass

    def _open_embedded_3d(self) -> None:
        pdb = (self.pdb_id.get() if hasattr(self, "pdb_id") else "").strip().upper()

        # If user didn't provide a PDB ID, try to use the selected RCSB hit (if any).
        if not pdb and hasattr(self, "vis_hits_list"):
            try:
                picked = self._vis_get_selected_pdb_id()
            except Exception:
                picked = ""
            if picked:
                pdb = picked
                try:
                    self.pdb_id.delete(0, "end")
                    self.pdb_id.insert(0, pdb)
                except Exception:
                    pass

        if not pdb:
            messagebox.showinfo(
                "No structure selected",
                "Mol* cannot visualize a raw DNA/protein sequence by itself — it needs a structure (PDB/mmCIF).\n\n"
                "Use the panel: 'Find a structure automatically (sequence similarity search)' to get a PDB ID,\n"
                "then select a hit and open it in Mol*.",
            )
            return

        if importlib.util.find_spec("webview") is None:
            messagebox.showinfo(
                "Embedded viewer not installed",
                "To use the embedded 3D viewer, install pywebview in the SAME Python that runs this app.\n\n"
                "If you are using the project virtual environment on Windows:\n"
                "  .\\.venv\\Scripts\\python.exe -m pip install pywebview\n\n"
                "Otherwise (global Python):\n"
                "  pip install pywebview\n\n"
                "Opening the browser viewer instead.",
            )
            url = f"https://molstar.org/viewer/?pdb={pdb}"
            webbrowser.open(url)
            return
        viewer_py = Path("viewer.py")
        if not viewer_py.exists():
            messagebox.showerror("Viewer missing", "viewer.py not found in project root.")
            return

        args = [sys.executable, str(viewer_py), "--pdb", pdb]

        try:
            subprocess.Popen(args, cwd=str(Path.cwd()))
        except Exception as e:
            messagebox.showerror(
                "Failed to open embedded viewer",
                f"Could not start embedded viewer.\n\nError: {e}\n\nFalling back to browser.",
            )
            url = f"https://molstar.org/viewer/?pdb={pdb}" if pdb else "https://molstar.org/viewer/"
            webbrowser.open(url)

    def _save_last_protein_fasta(self) -> None:
        protein = self.last_protein_alt or self.last_protein_ref or self.last_protein_single
        if not protein:
            messagebox.showinfo(
                "No protein available",
                "Run Single Analysis or Compare first so the app can extract a protein from an ORF.",
            )
            return

        path = filedialog.asksaveasfilename(
            title="Save protein FASTA",
            defaultextension=".fasta",
            filetypes=[("FASTA", "*.fasta"), ("All files", "*.*")],
        )
        if not path:
            return

        header = ">protein_from_app"
        wrapped = _format_wrapped(protein, 80)
        try:
            Path(path).write_text(header + "\n" + wrapped + "\n", encoding="utf-8")
        except Exception as e:
            messagebox.showerror("Save failed", f"Could not save file.\n\nError: {e}")
            return

        messagebox.showinfo("Saved", f"Saved: {path}")

    # ----------------------
    # Helpers
    # ----------------------
    def _set_text(self, widget: ScrolledText, text: str) -> None:
        widget.config(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", text)
        widget.config(state="disabled")


def main() -> None:
    app = DNAApp()

    # preload an example if present
    example = Path("examples/single_example.fasta")
    if example.exists():
        try:
            rec = read_fasta_file(example)
            app.single_input.insert("1.0", f">{rec.header}\n" if rec.header else "")
            app.single_input.insert("end", _format_wrapped(rec.sequence, 80))
            app.single_header.config(text=f"Header: {rec.header}" if rec.header else "")
        except Exception:
            pass

    app.mainloop()


if __name__ == "__main__":
    main()
