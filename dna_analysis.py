"""Core DNA sequence analysis functions."""
from collections import Counter
from math import log2

VALID_BASES = set("ACGT")


def normalize_sequence(sequence: str) -> str:
    return "".join(sequence.upper().split())


def is_valid_dna(sequence: str) -> bool:
    sequence = normalize_sequence(sequence)
    return bool(sequence) and set(sequence).issubset(VALID_BASES)


def nucleotide_counts(sequence: str) -> dict[str, int]:
    sequence = normalize_sequence(sequence)
    if not is_valid_dna(sequence):
        raise ValueError("Sequence must contain only A, C, G, and T.")
    counts = Counter(sequence)
    return {base: counts.get(base, 0) for base in "ACGT"}


def nucleotide_percentages(sequence: str) -> dict[str, float]:
    counts = nucleotide_counts(sequence)
    length = len(normalize_sequence(sequence))
    return {base: (count / length) * 100 for base, count in counts.items()}


def gc_content(sequence: str) -> float:
    counts = nucleotide_counts(sequence)
    length = len(normalize_sequence(sequence))
    return ((counts["G"] + counts["C"]) / length) * 100


def kmer_counts(sequence: str, k: int) -> dict[str, int]:
    sequence = normalize_sequence(sequence)
    if not is_valid_dna(sequence):
        raise ValueError("Sequence must contain only A, C, G, and T.")
    if k < 1 or k > len(sequence):
        raise ValueError("k must be between 1 and the sequence length.")
    return dict(Counter(sequence[i:i + k] for i in range(len(sequence) - k + 1)))


def shannon_entropy(sequence: str) -> float:
    counts = nucleotide_counts(sequence)
    length = len(normalize_sequence(sequence))
    entropy = 0.0
    for count in counts.values():
        if count:
            p = count / length
            entropy -= p * log2(p)
    return entropy


def sequence_summary(sequence: str) -> dict:
    sequence = normalize_sequence(sequence)
    counts = nucleotide_counts(sequence)
    result = {
        "length": len(sequence),
        "A_count": counts["A"],
        "C_count": counts["C"],
        "G_count": counts["G"],
        "T_count": counts["T"],
        "GC_content_percent": round(gc_content(sequence), 4),
        "Shannon_entropy_bits": round(shannon_entropy(sequence), 4),
    }
    for k in (1, 2, 3):
        kmers = kmer_counts(sequence, k)
        top_kmer, top_count = Counter(kmers).most_common(1)[0]
        result[f"unique_{k}mer_count"] = len(kmers)
        result[f"top_{k}mer"] = top_kmer
        result[f"top_{k}mer_count"] = top_count
    return result
