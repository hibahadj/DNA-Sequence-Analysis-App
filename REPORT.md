# Mini‑Project Report — DNA Sequence Analysis & Mutation Detection (GUI)

Date: February 2026

## 0) Executive summary
This project delivers a small bioinformatics application with a graphical interface to:

1) Analyze a single DNA sequence (validation, base statistics, reverse‑complement, codons, ORFs, translation).
2) Compare two DNA sequences (global alignment, mutation detection + classification, mutation rate, protein/ORF impact summary).
3) Help the user visualize and interpret results (colored alignment, mutation table, CSV export, copyable report) and connect the output to common online 3D viewers (Mol*, RCSB PDB, AlphaFold DB).

The application is written in Python and uses Tkinter (standard library) for a lightweight Windows‑friendly GUI.

---

## 1) What the assignment requires (mapping to delivered features)
Based on the provided TP exercises (TP1→TP4) and the mini‑project description, the typical requirements are:

- Read DNA from FASTA or direct input.
- Validate DNA alphabet (A/T/C/G).
- Compute basic statistics and transformations.
- Search/scan biological features (codons, ORFs, translations).
- Compare two sequences, align them, and detect mutations.
- Produce an interpretable output (summary + visualization).
- Provide a clear report explaining the theory, the workflow, and the tools/functions used.

This application implements all of the above in a single integrated GUI.

---

## 2) Tools and technologies used

### 2.1 Languages, environment
- Python 3.x
- OS target: Windows (but code is cross‑platform)

### 2.2 Libraries
- Tkinter / ttk / ScrolledText (standard library): GUI
- textwrap, pathlib, webbrowser, subprocess, importlib.util, csv, random (standard library): utilities
- Optional: `pywebview` for embedded 3D viewer window

### 2.3 Online tools (used for visualization)
- Mol* Viewer (web): interactive 3D molecular viewer
- RCSB PDB: protein structure database
- AlphaFold DB: predicted structures (when available)

---

## 3) User guide — how to use the app (step‑by‑step)

### 3.1 Run
From the project root:

```bash
python app.py
```

### 3.2 Inputs supported
In all text areas you can paste:
- FASTA text (optional header line starting with `>`)
- Raw DNA (A/T/C/G) directly

Spaces/newlines are ignored during normalization; only A/T/C/G are accepted.

### 3.3 Single Analysis tab
Goal: analyze 1 DNA sequence.

Steps:
1) Click “Load FASTA file” or “Load example” (or paste your sequence).
2) Optional: enable “Also scan reverse‑complement (6‑frame ORFs)”.
3) Click “Analyze”.

Output includes:
- Length and nucleotide counts
- GC% and AT%
- Reverse‑complement
- Codon statistics (start/stop counts)
- ORF list (first 10 shown)
- Protein translation for 3 reading frames

### 3.4 Compare tab
Goal: align REF vs MUT and detect mutations.

Steps:
1) Load two sequences using “Load REF FASTA” and “Load MUT FASTA”.
	 - Or click “Load random example” to quickly test.
2) Click “Compare”.

What you get:
- Mutation summary (counts of substitutions/insertions/deletions)
- Transition/transversion breakdown
- Mutation rate per reference base
- ORF comparison (+ strand)
- Protein impact (heuristic, longest ORF)
- Alignment visualization
	- Bases are colored (match/transition/transversion/insertion/deletion)
	- A small legend is displayed above the output

Extra productivity tools:
- “Mutation table”: opens a table (type, positions, bases, class)
- “Export CSV”: saves the mutation list to a CSV file
- “Copy report”: copies the entire Compare output text to the clipboard

### 3.5 Visualization tab (3D)
Goal: open online tools for protein 3D visualization.

Two main use cases:

**A) You already have a PDB ID**
1) Go to the Visualization tab.
2) Enter the PDB ID (example: `1CRN`).
3) Click “Open embedded 3D”.

If embedded view does not work, the app automatically falls back to the browser.

**B) You do not have a PDB ID**
1) Run Single Analysis or Compare first (so the app can extract a protein from an ORF).
2) In Visualization, click “Save last protein FASTA”.
3) Use the saved protein FASTA to search databases/tools outside the app.

### 3.6 Guide tab
The Guide tab is an in‑app manual:
- Explains the alignment coloring and what each mutation means
- Provides step‑by‑step instructions for Mol*/RCSB/AlphaFold usage

---

## 4) Theory (mutations + alignment + translation)

### 4.1 DNA basics
DNA is composed of 4 nucleotides: A, T, C, G.

Complementary base pairing:
- A pairs with T
- C pairs with G

Reverse‑complement is used because genes can be encoded on either strand.

### 4.2 Codons, reading frames, and translation
Translation reads DNA in groups of three bases (codons). Since codons are triplets, there are 3 reading frames on a strand:

- Frame 0: start at position 1
- Frame 1: start at position 2
- Frame 2: start at position 3

Stop codons (in DNA): TAA, TAG, TGA.
The app translates codons using the standard genetic code (unknown codon → `X`).

### 4.3 ORFs (Open Reading Frames)
An ORF is a region that can be translated into protein. In this project’s simplified educational definition:

- Start codon: ATG
- ORF ends at the first in‑frame stop codon

ORFs are scanned in each frame; optionally the reverse‑complement is also scanned.

### 4.4 Mutation types (DNA level)
When comparing a reference (REF) to a mutated/alternate sequence (MUT):

- Substitution: one base replaced by another (e.g., A→G)
- Insertion: an extra base appears in MUT relative to REF
- Deletion: a base is missing in MUT relative to REF

Substitutions are further classified into:

- Transition: A↔G or C↔T
- Transversion: any other substitution (A↔C, A↔T, G↔C, G↔T)

Why transitions/transversions matter: transitions are often more frequent biologically because they substitute within the same nucleotide class (purine↔purine or pyrimidine↔pyrimidine).

### 4.5 Frameshift, missense, nonsense (protein‑level impact)
If an insertion/deletion changes the coding sequence length by a value not divisible by 3, the reading frame changes.
This is a **frameshift**, and it usually changes many amino acids downstream.

Substitutions can be:
- Silent (synonymous): amino acid does not change
- Missense: amino acid changes
- Nonsense: substitution creates a stop codon, shortening the protein

In this project, protein impact is estimated with a best‑effort heuristic by comparing the longest ORF translated proteins.

### 4.6 Why we need alignment
Two sequences may differ by insertions/deletions, which shift indices. Comparing raw strings position‑by‑position is not reliable.
Alignment introduces gaps to maximize similarity and place mutations in the correct context.

### 4.7 Global alignment (Needleman–Wunsch) — dynamic programming
This project uses a global alignment algorithm suitable for aligning full‑length sequences.
With scoring parameters:

- match = +1
- mismatch = −1
- gap = −1

Let $dp[i][j]$ be the best score aligning the first $i$ characters of REF and the first $j$ characters of MUT.

Recurrence:
$$
dp[i][j] = \max\begin{cases}
dp[i-1][j-1] + s(\text{ref}[i], \text{alt}[j])\\
dp[i-1][j] + gap\\
dp[i][j-1] + gap
\end{cases}
$$
Then a traceback reconstructs the aligned strings.

Time complexity: $O(nm)$, space complexity: $O(nm)$.

### 4.8 Mutation rate
Mutation rate reported is:
$$\text{rate} = \frac{\#\text{mutations}}{\text{length of REF}}$$

This is a simple per‑base rate for an educational comparison.

---

## 5) Workflow (steps of work)
1) Implement reusable bioinformatics functions in a small package (`dna_app/`).
2) Add FASTA parsing for both file input and pasted text.
3) Implement global alignment (Needleman–Wunsch) and mutation extraction from the alignment.
4) Build Tkinter GUI:
	 - Single analysis
	 - Compare + visualization
	 - Visualization links
5) Improve UI/UX:
	 - Colored mutations in the alignment
	 - Mutation table + export/copy
	 - Guide tab explaining how to read results and use online tools

---

## 6) Implementation details (modules and main functions)

### 6.1 Core package `dna_app/`

**FASTA** (`dna_app/fasta.py`)
- `parse_fasta_text(text) -> FastaRecord`: accepts FASTA or plain DNA
- `read_fasta_file(path) -> FastaRecord`: reads FASTA file from disk

**DNA utilities** (`dna_app/dna.py`)
- `validate_dna(seq)`: checks only A/T/C/G and non‑empty
- `stats(seq) -> SequenceStats`: length, A/T/C/G counts, GC%, AT%
- `reverse_complement(seq)`: reverse‑complement transformation
- `split_codons(seq, frame)`: codon list for a reading frame
- `translate_dna(seq, frame, stop_at_stop)`: translate into amino acids
- `find_orfs(seq, include_reverse_complement) -> list[ORF]`: ORF detection per frame

**Alignment** (`dna_app/alignment.py`)
- `needleman_wunsch_global(ref, alt, match=1, mismatch=-1, gap=-1) -> AlignmentResult`:
	returns `aligned_ref`, `aligned_alt`, `midline`, `score`

**Mutations** (`dna_app/mutations.py`)
- `detect_mutations_from_alignment(aln) -> list[Mutation]`: generates substitutions/insertions/deletions with 1‑based coordinates
- `classify_substitution(ref_base, alt_base)`: transition/transversion/other
- `summarize_mutations(mutations)`: counts per class
- `mutation_rate(mutations, ref_length)`: per‑base rate

**Protein impact (heuristic)** (`dna_app/impact.py`)
- `analyze_protein_impact(ref_dna, alt_dna) -> ProteinImpact`:
	compares longest ORF proteins to label likely `silent/missense/nonsense/frameshift/unknown`

### 6.2 GUI application

**Main GUI** (`app.py`)
- Tabs:
	- Single Analysis
	- Compare
	- Visualization
	- Guide
- Key GUI methods (selected):
	- `_single_analyze()`: runs `validate_dna`, `stats`, `reverse_complement`, `find_orfs`, `translate_dna`
	- `_compare_run()`: runs `needleman_wunsch_global`, `detect_mutations_from_alignment`, `summarize_mutations`, `mutation_rate`, `analyze_protein_impact`
	- `_show_mutation_table()`: displays mutations in a table
	- `_export_mutations_csv()`: exports mutations

**Embedded 3D window** (`viewer.py`)
- Uses `pywebview` (optional) to embed Mol* viewer.

---

## 7) How to reproduce results (examples)

### 7.1 Quick demo
1) Run the app.
2) Go to Compare.
3) Click “Load random example”.
4) Click “Compare”.
5) Open “Mutation table” and export CSV.

### 7.2 Provided example files
- `examples/single_example.fasta`
- `examples/reference.fasta` and `examples/mutated.fasta`
- `TP4/origin.fasta` and `TP4/mutation.fasta`

---

## 8) Limitations and future improvements
- Protein impact is a heuristic (longest ORF only); a full solution would map each DNA mutation to codon/residue positions.
- Global alignment is $O(nm)$; very long sequences may be slow.
- 3D visualization is external (Mol*/RCSB/AlphaFold). Automatic “highlight mutated residues in 3D” requires mapping DNA→protein indices and a structure.

Future improvements:
- Add residue‑level mapping from mutations to protein positions (for ORF region)
- Add export of the colored alignment as HTML (easy to submit online)
- Add local alignment (Smith–Waterman) as an option for partial matching

---

## 9) Conclusion
The delivered application unifies the main DNA analysis and mutation detection tasks into one clean GUI workflow, with theory‑aligned outputs, clear mutation visualization, and practical export options suitable for a course mini‑project submission.
