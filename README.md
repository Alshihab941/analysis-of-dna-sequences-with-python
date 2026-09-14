# Analysis of DNA Sequences with Python

This repository contains the computational foundation for the research project **Analysis of DNA Sequences with Python**.

## Current computational foundation
- DNA sequence representation and validation
- Nucleotide composition (A, C, G, T)
- GC content
- k-mer analysis (k=1, 2, 3)
- Shannon entropy as a sequence-complexity measure
- FASTA input using Biopython `Bio.SeqIO`
- CSV output for downstream analysis

## Structure
```text
analysis-of-dna-sequences-with-python/
├── README.md
├── requirements.txt
├── data/
│   └── sequences.fasta
├── src/
│   ├── dna_analysis.py
│   └── analyze_fasta.py
├── results/
│   └── analysis_results.csv
├── references/
│   └── resources.md
└── research_stage_notes.md
```

## Setup
Use Python 3.10+.

```bash
python -m pip install -r requirements.txt
```

## Run
From the repository root:

```bash
python src/analyze_fasta.py
```

The program reads `data/sequences.fasta`, analyzes each record, prints a summary, and writes `results/analysis_results.csv`.

## Methods
### Nucleotide composition
Counts A, C, G, and T and reports their percentages.

### GC content
`(G + C) / sequence_length * 100`

### k-mer analysis
A k-mer is a substring of length k. The pipeline counts k-mers for k=1, 2, and 3 and reports the most frequent one.

### Shannon entropy
`H = -Σ p_i log2(p_i)`

The calculation measures diversity in the nucleotide probability distribution.

## Research-stage note
The included FASTA sequences are **synthetic demonstration sequences**. They are only for testing the computational workflow and are not biological findings. The final biological dataset, research question, and statistical analysis should be selected after discussion with the professor.
