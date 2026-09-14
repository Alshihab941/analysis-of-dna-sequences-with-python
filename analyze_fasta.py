"""Reproducible FASTA analysis pipeline."""
from pathlib import Path
import csv
from Bio import SeqIO
from dna_analysis import sequence_summary

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_ROOT / "data" / "sequences.fasta"
OUTPUT_FILE = PROJECT_ROOT / "results" / "analysis_results.csv"


def analyze_fasta(input_file: Path):
    rows = []
    for record in SeqIO.parse(input_file, "fasta"):
        metrics = sequence_summary(str(record.seq))
        metrics["sequence_id"] = record.id
        metrics["description"] = record.description
        rows.append(metrics)
    return rows


def write_csv(rows, output_file: Path):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with output_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main():
    rows = analyze_fasta(INPUT_FILE)
    write_csv(rows, OUTPUT_FILE)
    print("DNA sequence analysis completed.")
    print(f"Input:  {INPUT_FILE}")
    print(f"Output: {OUTPUT_FILE}")
    for row in rows:
        print(f"{row['sequence_id']}: length={row['length']}, GC={row['GC_content_percent']:.2f}%, entropy={row['Shannon_entropy_bits']:.4f} bits")


if __name__ == "__main__":
    main()
