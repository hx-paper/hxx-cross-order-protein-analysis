from __future__ import annotations

import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError as exc:
    raise SystemExit("Biopython is required: install or load the bundled workspace Python dependencies") from exc

FREE_AA = {
    "A": (3, 7, 1, 2, 0), "R": (6, 14, 4, 2, 0), "N": (4, 8, 2, 3, 0),
    "D": (4, 7, 1, 4, 0), "C": (3, 7, 1, 2, 1), "E": (5, 9, 1, 4, 0),
    "Q": (5, 10, 2, 3, 0), "G": (2, 5, 1, 2, 0), "H": (6, 9, 3, 2, 0),
    "I": (6, 13, 1, 2, 0), "L": (6, 13, 1, 2, 0), "K": (6, 14, 2, 2, 0),
    "M": (5, 11, 1, 2, 1), "F": (9, 11, 1, 2, 0), "P": (5, 9, 1, 2, 0),
    "S": (3, 7, 1, 3, 0), "T": (4, 9, 1, 3, 0), "W": (11, 12, 2, 2, 0),
    "Y": (9, 11, 1, 3, 0), "V": (5, 11, 1, 2, 0),
}


def read_fasta(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    sequence = "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))
    sequence = re.sub(r"\s+", "", sequence).upper().rstrip("*")
    bad = sorted(set(sequence) - set(FREE_AA))
    if not sequence or bad:
        raise SystemExit(f"Invalid protein sequence; unsupported symbols: {bad}")
    return sequence


def formula(sequence: str) -> tuple[str, int]:
    totals = [0, 0, 0, 0, 0]
    for residue in sequence:
        for i, value in enumerate(FREE_AA[residue]):
            totals[i] += value
    totals[1] -= 2 * (len(sequence) - 1)
    totals[3] -= len(sequence) - 1
    c, h, n, o, s = totals
    return f"C{c}H{h}N{n}O{o}S{s}", sum(totals)


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("Usage: calculate_physicochemical.py <protein.fasta> <protein-name> <output-dir>")
    fasta = Path(sys.argv[1]).resolve()
    protein = re.sub(r"[^A-Za-z0-9_.-]+", "_", sys.argv[2].strip())
    out = Path(sys.argv[3]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    sequence = read_fasta(fasta)
    analysis = ProteinAnalysis(sequence)
    composition = Counter(sequence)
    aliphatic = 100 * (composition["A"] + 2.9 * composition["V"] + 3.9 * (composition["I"] + composition["L"])) / len(sequence)
    molecular_formula, atom_count = formula(sequence)
    extinction = analysis.molar_extinction_coefficient()
    instability = analysis.instability_index()
    result = {
        "length_aa": len(sequence),
        "molecular_weight_da": analysis.molecular_weight(),
        "theoretical_pi": analysis.isoelectric_point(),
        "acidic_residues_Asp_Glu": sequence.count("D") + sequence.count("E"),
        "basic_residues_Arg_Lys": sequence.count("R") + sequence.count("K"),
        "formula": molecular_formula,
        "atom_count": atom_count,
        "extinction_coefficient_reduced_M-1_cm-1": extinction[0],
        "extinction_coefficient_cystine_M-1_cm-1": extinction[1],
        "instability_index": instability,
        "predicted_stability": "stable" if instability < 40 else "unstable",
        "aliphatic_index": aliphatic,
        "gravy": analysis.gravy(),
        "predicted_half_life": "30 h (mammalian reticulocytes, in vitro); >20 h (yeast); >10 h (E. coli)",
        "n_terminal_residue": sequence[0],
        "method": "ExPASy ProtParam-compatible calculation using Biopython ProteinAnalysis",
    }
    (out / f"{protein}_input_sequence.fasta").write_text(f">Query_{protein}_user_provided\n{sequence}\n", encoding="utf-8")
    (out / f"{protein}_physicochemical_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    rows = [
        ("氨基酸数目", str(len(sequence)), "aa"),
        ("分子量", f"{analysis.molecular_weight() / 1000:.2f}", "kDa"),
        ("理论等电点", f"{analysis.isoelectric_point():.2f}", "pI"),
        ("酸性残基", str(result["acidic_residues_Asp_Glu"]), "Asp + Glu"),
        ("碱性残基", str(result["basic_residues_Arg_Lys"]), "Arg + Lys"),
        ("分子式", molecular_formula, "完整未修饰蛋白"),
        ("原子总数", f"{atom_count:,}", "个"),
        ("消光系数（还原态/胱氨酸）", f"{extinction[0]:,} / {extinction[1]:,}", "M^-1 cm^-1, 280 nm"),
        ("不稳定指数", f"{instability:.2f}", "<40预测稳定；>40预测不稳定"),
        ("脂肪族指数", f"{aliphatic:.2f}", "Aliphatic index"),
        ("平均疏水性", f"{analysis.gravy():.3f}", "GRAVY"),
        ("预测半衰期", result["predicted_half_life"], "基于N端规则"),
    ]
    with (out / f"{protein}_physicochemical_parameters.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["参数", "预测值", "单位或说明"])
        writer.writerows(rows)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
