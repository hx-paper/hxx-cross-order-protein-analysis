from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from PIL import Image
from scipy import ndimage

from convert_modelcif_to_pdb import parse_atom_site


CLASS_NAMES = ["disallowed", "generously_allowed", "additionally_allowed", "most_favoured"]
CLASSIC = ["#FFFFFF", "#FFF9A6", "#FFF200", "#F50000"]
PUBLICATION = ["#FAFAFA", "#F4E6A2", "#B9DCD5", "#5FA8A0"]


def cif_atoms(path: Path) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    for atom in parse_atom_site(path):
        if atom["group_PDB"] != "ATOM":
            continue
        normalized.append({
            "chain": atom.get("auth_asym_id", atom.get("label_asym_id", "A")),
            "seq": int(atom.get("auth_seq_id", atom["label_seq_id"])),
            "ins": atom.get("pdbx_PDB_ins_code", "."),
            "resname": atom["label_comp_id"],
            "atom": atom["label_atom_id"],
            "xyz": np.array([float(atom["Cartn_x"]), float(atom["Cartn_y"]), float(atom["Cartn_z"])]),
            "bfactor": float(atom.get("B_iso_or_equiv", 0.0)),
        })
    return normalized


def pdb_atoms(path: Path) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if not line.startswith("ATOM  "):
            continue
        normalized.append({
            "chain": line[21:22].strip() or "A",
            "seq": int(line[22:26]),
            "ins": line[26:27].strip() or ".",
            "resname": line[17:20].strip(),
            "atom": line[12:16].strip(),
            "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            "bfactor": float(line[60:66] or 0.0),
        })
    return normalized


def load_atoms(path: Path) -> list[dict[str, object]]:
    atoms = pdb_atoms(path) if path.suffix.lower() in {".pdb", ".ent"} else cif_atoms(path)
    if not atoms:
        raise ValueError(f"No ATOM coordinates found in {path}")
    return atoms


def dihedral(p0, p1, p2, p3) -> float:
    p0, p1, p2, p3 = map(lambda point: np.asarray(point, dtype=float), (p0, p1, p2, p3))
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def compute_angles(atoms: list[dict[str, object]]) -> list[dict[str, object]]:
    residues: dict[tuple[str, int, str, str], dict[str, object]] = defaultdict(lambda: {"atoms": {}, "bfactors": []})
    for atom in atoms:
        key = (str(atom["chain"]), int(atom["seq"]), str(atom["ins"]), str(atom["resname"]))
        residues[key]["atoms"][str(atom["atom"])] = atom["xyz"]
        residues[key]["bfactors"].append(float(atom["bfactor"]))
    chains: dict[str, list[dict[str, object]]] = defaultdict(list)
    for (chain, seq, ins, name), data in residues.items():
        chains[chain].append({"chain": chain, "seq": seq, "ins": ins, "name": name, **data})
    for chain in chains:
        chains[chain].sort(key=lambda residue: (residue["seq"], residue["ins"]))

    rows: list[dict[str, object]] = []
    for chain, chain_residues in sorted(chains.items()):
        for index in range(1, len(chain_residues) - 1):
            previous, current, following = chain_residues[index - 1:index + 2]
            required = [
                previous["atoms"].get("C"), current["atoms"].get("N"), current["atoms"].get("CA"),
                current["atoms"].get("C"), following["atoms"].get("N"),
            ]
            if any(value is None for value in required):
                continue
            if np.linalg.norm(previous["atoms"]["C"] - current["atoms"]["N"]) > 2.0:
                continue
            if np.linalg.norm(current["atoms"]["C"] - following["atoms"]["N"]) > 2.0:
                continue
            phi = dihedral(previous["atoms"]["C"], current["atoms"]["N"], current["atoms"]["CA"], current["atoms"]["C"])
            psi = dihedral(current["atoms"]["N"], current["atoms"]["CA"], current["atoms"]["C"], following["atoms"]["N"])
            residue_class = "glycine" if current["name"] == "GLY" else "proline" if current["name"] == "PRO" else "general"
            rows.append({
                "chain": chain,
                "residue_number": int(current["seq"]),
                "insertion_code": "" if current["ins"] in {"?", "."} else current["ins"],
                "residue_name": current["name"],
                "residue_class": residue_class,
                "phi": phi,
                "psi": psi,
                "confidence": float(np.mean(current["bfactors"])),
            })
    if not rows:
        raise ValueError("No valid phi/psi pairs could be calculated")
    return rows


def official_region_map(path: Path) -> np.ndarray:
    rgb = np.asarray(Image.open(path).convert("RGB"))
    height, width = rgb.shape[:2]
    if height < 800 or width < 600:
        raise ValueError("Official PROCHECK PNG is too small for reliable region recovery")
    # SAVES/PROCHECK main plot occupies these stable fractions of the portrait output.
    y0, y1 = round(height * 0.158), round(height * 0.6285)
    x0, x1 = round(width * 0.1645), round(width * 0.8445)
    crop = rgb[y0:y1, x0:x1]
    if abs(crop.shape[0] - crop.shape[1]) / max(crop.shape[:2]) > 0.08:
        raise ValueError(f"Automatic PROCHECK plot crop is not square enough: {crop.shape[:2]}")
    red, green, blue = crop[..., 0], crop[..., 1], crop[..., 2]
    region = np.full(crop.shape[:2], -1, dtype=np.int8)
    region[(red > 225) & (green > 225) & (blue > 225)] = 0
    region[(red > 180) & (green > 180) & (blue >= 80) & (blue < 225)] = 1
    region[(red > 180) & (green > 180) & (blue < 80)] = 2
    region[(red > 150) & (green < 100) & (blue < 100)] = 3
    minimum_component = max(100, round(region.size * 0.0005))
    for category in (1, 2, 3):
        labelled, _ = ndimage.label(region == category)
        sizes = np.bincount(labelled.ravel())
        keep = sizes >= minimum_component
        keep[0] = False
        region[(region == category) & ~keep[labelled]] = -1
    invalid = region < 0
    if invalid.all():
        raise ValueError("Could not recover colored PROCHECK regions from the official PNG")
    nearest = ndimage.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    region[invalid] = region[tuple(nearest[:, invalid])]
    yi = np.linspace(0, region.shape[0] - 1, 361).round().astype(int)
    xi = np.linspace(0, region.shape[1] - 1, 361).round().astype(int)
    return np.flipud(region[np.ix_(yi, xi)])


def parse_summary(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        r"Ramachandran plot:\s*([\d.]+)%\s+core\s+([\d.]+)%\s+allow\s+([\d.]+)%\s+gener\s+([\d.]+)%\s+disall",
        text,
        flags=re.IGNORECASE,
    )
    if not match:
        raise ValueError("Could not parse the official PROCHECK Ramachandran summary line")
    job_match = re.search(r"/Jobs/(\d+)/", text, flags=re.IGNORECASE)
    residues_match = re.search(r"\b(\d+)\s+residues\s*\|", text, flags=re.IGNORECASE)
    return {
        "percentages": {
            "most_favoured": float(match.group(1)),
            "additionally_allowed": float(match.group(2)),
            "generously_allowed": float(match.group(3)),
            "disallowed": float(match.group(4)),
        },
        "saves_job": job_match.group(1) if job_match else "not recorded",
        "total_residues": int(residues_match.group(1)) if residues_match else None,
    }


def assign_regions(rows: list[dict[str, object]], region_map: np.ndarray) -> None:
    for row in rows:
        x = int(np.clip(round(float(row["phi"]) + 180), 0, 360))
        y = int(np.clip(round(float(row["psi"]) + 180), 0, 360))
        category = int(region_map[y, x])
        row["official_region"] = CLASS_NAMES[category]


def official_counts(percentages: dict[str, float], denominator: int) -> dict[str, int]:
    counts = {
        "additionally_allowed": round(denominator * percentages["additionally_allowed"] / 100),
        "generously_allowed": round(denominator * percentages["generously_allowed"] / 100),
        "disallowed": round(denominator * percentages["disallowed"] / 100),
    }
    counts["most_favoured"] = denominator - sum(counts.values())
    return counts


def residue_label(row: dict[str, object]) -> str:
    return f"{row['chain']}:{str(row['residue_name']).title()}{row['residue_number']}"


def draw_scatter(ax, rows: list[dict[str, object]], style: str) -> None:
    if style == "classic":
        specifications = {
            "general": ("s", 15, "#000000", "#000000"),
            "glycine": ("^", 25, "#000000", "#000000"),
            "proline": ("s", 23, "#7A327A", "none"),
        }
    else:
        specifications = {
            "general": ("o", 17, "#344054", "#FFFFFF"),
            "glycine": ("^", 27, "#246BCE", "#FFFFFF"),
            "proline": ("s", 23, "#7A5195", "#FFFFFF"),
        }
    for residue_class, (marker, size, edge, face) in specifications.items():
        subset = [row for row in rows if row["residue_class"] == residue_class]
        ax.scatter(
            [row["phi"] for row in subset], [row["psi"] for row in subset], marker=marker,
            s=size, facecolor=face, edgecolor=edge, linewidth=0.65, alpha=0.92, zorder=4,
        )
    generous = [row for row in rows if row["residue_class"] == "general" and row["official_region"] == "generously_allowed"]
    if generous:
        ax.scatter([row["phi"] for row in generous], [row["psi"] for row in generous], s=54,
                   facecolor="none", edgecolor="#D97706", linewidth=1.1, zorder=6)
    disallowed = [row for row in rows if row["residue_class"] == "general" and row["official_region"] == "disallowed"]
    for index, row in enumerate(disallowed[:8]):
        ax.scatter(row["phi"], row["psi"], s=68, marker="D", facecolor="#D92D20",
                   edgecolor="white", linewidth=0.8, zorder=7)
        dx = 14 if float(row["phi"]) < 130 else -14
        horizontal = "left" if dx > 0 else "right"
        ax.annotate(
            residue_label(row), (row["phi"], row["psi"]), xytext=(dx, 10 + 10 * (index % 2)),
            textcoords="offset points", fontsize=8, color="#B42318", ha=horizontal, va="center",
            arrowprops={"arrowstyle": "-", "color": "#B42318", "lw": 0.75}, zorder=8,
        )


def render_plot(
    region_map: np.ndarray,
    rows: list[dict[str, object]],
    percentages: dict[str, float],
    counts: dict[str, int],
    denominator: int,
    palette: list[str],
    style: str,
    title: str,
    subtitle: str,
    output_stem: Path,
) -> None:
    fig = plt.figure(figsize=(7.4, 7.6), facecolor="white")
    grid = fig.add_gridspec(2, 1, height_ratios=[12, 1.25], hspace=0.31,
                            left=0.13, right=0.97, top=0.88, bottom=0.16)
    ax = fig.add_subplot(grid[0, 0])
    bar = fig.add_subplot(grid[1, 0])
    ax.imshow(region_map, origin="lower", extent=(-180, 180, -180, 180),
              cmap=ListedColormap(palette), interpolation="nearest", vmin=-0.5, vmax=3.5, zorder=0)
    ax.contour(np.linspace(-180, 180, 361), np.linspace(-180, 180, 361), region_map,
               levels=[0.5, 1.5, 2.5], colors="#52606D", linewidths=0.45, alpha=0.72, zorder=1)
    ax.axhline(0, color="#667085", lw=0.45, alpha=0.55, zorder=2)
    ax.axvline(0, color="#667085", lw=0.45, alpha=0.55, zorder=2)
    draw_scatter(ax, rows, style)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-180, 181, 60))
    ax.set_xlabel(r"$\phi$ (°)", fontsize=10)
    ax.set_ylabel(r"$\psi$ (°)", fontsize=10)
    ax.tick_params(direction="out", length=3, width=0.7, labelsize=9)
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color("#344054")
    ax.set_title(f"{title}\n{subtitle}", loc="left", fontsize=12.5, fontweight="bold", pad=9)

    ordered_keys = ["most_favoured", "additionally_allowed", "generously_allowed", "disallowed"]
    ordered_colors = [palette[3], palette[2], palette[1], palette[0]]
    left = 0.0
    for key, color in zip(ordered_keys, ordered_colors):
        value = percentages[key]
        bar.barh(0, value, left=left, height=0.52, color=color, edgecolor="white", linewidth=0.8)
        if value >= 7:
            text_color = "white" if key == "most_favoured" else "#24313A"
            bar.text(left + value / 2, 0, f"{value:.1f}%", ha="center", va="center",
                     fontsize=8.2, color=text_color, fontweight="bold")
        left += value
    bar.set_xlim(0, 100)
    bar.set_ylim(-0.65, 0.65)
    bar.set_yticks([])
    bar.set_xticks([0, 25, 50, 75, 100], ["0", "25", "50", "75", "100%"])
    bar.tick_params(axis="x", labelsize=8, length=2.5, width=0.6, pad=2)
    bar.set_title(
        f"Official PROCHECK regions (non-Gly/Pro, n={denominator}): "
        f"{counts['most_favoured']} / {counts['additionally_allowed']} / "
        f"{counts['generously_allowed']} / {counts['disallowed']}",
        fontsize=8.4, loc="left", pad=3, color="#344054",
    )
    for side in ("top", "left", "right"):
        bar.spines[side].set_visible(False)
    bar.spines["bottom"].set_color("#98A2B3")
    bar.spines["bottom"].set_linewidth(0.6)

    labels = ["Most favoured", "Additionally allowed", "Generously allowed", "Disallowed"]
    region_handles = [Patch(facecolor=color, edgecolor="#667085", linewidth=0.5, label=label)
                      for color, label in zip(ordered_colors, labels)]
    if style == "classic":
        marker_handles = [
            Line2D([], [], marker="s", linestyle="none", markerfacecolor="black", markeredgecolor="black", markersize=4.5, label="Other residue"),
            Line2D([], [], marker="^", linestyle="none", markerfacecolor="black", markeredgecolor="black", markersize=5.5, label="Gly"),
            Line2D([], [], marker="s", linestyle="none", markerfacecolor="none", markeredgecolor="#7A327A", markersize=5, label="Pro"),
        ]
    else:
        marker_handles = [
            Line2D([], [], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor="#344054", markersize=5, label="Other residue"),
            Line2D([], [], marker="^", linestyle="none", markerfacecolor="white", markeredgecolor="#246BCE", markersize=5.5, label="Gly"),
            Line2D([], [], marker="s", linestyle="none", markerfacecolor="white", markeredgecolor="#7A5195", markersize=5, label="Pro"),
        ]
    fig.legend(handles=region_handles + marker_handles, loc="lower center", bbox_to_anchor=(0.5, 0.045),
               ncol=4, frameon=False, fontsize=8.1, handlelength=1.3, columnspacing=1.1, handletextpad=0.4)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output_stem.with_suffix(f".{suffix}"), dpi=600 if suffix == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_thesis_text(
    output: Path,
    protein: str,
    stats: dict[str, object],
    counts: dict[str, int],
    denominator: int,
    disallowed: list[str],
    saves_version: str,
    procheck_version: str,
) -> None:
    p = stats["percentages"]
    disallowed_text = "、".join(disallowed) if disallowed else "未从主图坐标中检出可确认的禁阻区残基"
    text = f"""# {protein}蛋白Ramachandran图分析（简短版）

## 材料与方法

将预测蛋白结构转换为标准PDB格式，并提交至UCLA SAVES {saves_version}平台，采用PROCHECK {procheck_version}评价主链二面角分布。PROCHECK主图以非甘氨酸、非脯氨酸残基为统计对象。根据官方区域边界和统计结果，从结构坐标重新计算φ/ψ角，并分别绘制红—黄经典版和青绿—浅黄期刊版Ramachandran图；两图使用完全相同的数据，仅配色不同。

## 结果与分析

在{denominator}个参与PROCHECK主图统计的残基中，{counts['most_favoured']}个（{p['most_favoured']:.1f}%）位于最有利区，{counts['additionally_allowed']}个（{p['additionally_allowed']:.1f}%）位于额外允许区，{counts['generously_allowed']}个（{p['generously_allowed']:.1f}%）位于宽松允许区，{counts['disallowed']}个（{p['disallowed']:.1f}%）位于禁阻区。禁阻区残基为：{disallowed_text}。结果反映了该模型主链构象的立体化学合理性；对于禁阻区残基或低置信度区域，后续结构功能解释仍需谨慎，且本分析不能替代实验结构验证。

## 图注

**图D {protein}蛋白Ramachandran图（PROCHECK经典配色）。** 红色、黄色、浅黄色和白色区域分别表示最有利区、额外允许区、宽松允许区和禁阻区；方形、三角形及空心方形分别表示普通残基、甘氨酸和脯氨酸。下方堆叠条表示官方PROCHECK区域比例。

**图E {protein}蛋白Ramachandran图（期刊配色）。** 深青绿色、浅青绿色、浅黄色和白色区域分别表示最有利区、额外允许区、宽松允许区和禁阻区。该图与图D使用相同的φ/ψ坐标和官方统计，仅采用适合论文排版的替代配色。

> SAVES任务编号：{stats['saves_job']}。边界说明：Ramachandran图评价的是模型的立体化学合理性，不能单独证明预测结构即为蛋白的真实生物学构象。
"""
    output.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Render classic and publication Ramachandran plots from official PROCHECK outputs")
    parser.add_argument("structure", type=Path)
    parser.add_argument("official_plot", type=Path)
    parser.add_argument("summary", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("protein_name")
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--analysis-dir", type=Path)
    parser.add_argument("--saves-version", default="v6.1")
    parser.add_argument("--procheck-version", default="v3.5.4")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = args.data_dir or args.output_dir
    analysis_dir = args.analysis_dir or args.output_dir
    data_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir.mkdir(parents=True, exist_ok=True)
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", args.protein_name.strip())
    rows = compute_angles(load_atoms(args.structure.resolve()))
    region_map = official_region_map(args.official_plot.resolve())
    assign_regions(rows, region_map)
    stats = parse_summary(args.summary.resolve())
    general_rows = [row for row in rows if row["residue_class"] == "general"]
    denominator = len(general_rows)
    counts = official_counts(stats["percentages"], denominator)
    map_counts = {name: sum(row["official_region"] == name for row in general_rows) for name in CLASS_NAMES}
    discrepancies = {name: map_counts[name] - counts[name] for name in CLASS_NAMES}
    disallowed_rows = [row for row in general_rows if row["official_region"] == "disallowed"]
    disallowed_labels = [residue_label(row) for row in disallowed_rows]

    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 10, "axes.linewidth": 0.8,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })
    subtitle = f"SAVES job {stats['saves_job']} | {stats['total_residues'] or 'unknown'} residues"
    classic_stem = args.output_dir / f"Figure_D_{safe_name}_Ramachandran_classic_red_yellow_10pt"
    publication_stem = args.output_dir / f"Figure_E_{safe_name}_Ramachandran_publication_teal_yellow_10pt"
    render_plot(region_map, rows, stats["percentages"], counts, denominator, CLASSIC, "classic",
                f"{args.protein_name} Ramachandran plot", subtitle, classic_stem)
    render_plot(region_map, rows, stats["percentages"], counts, denominator, PUBLICATION, "publication",
                f"{args.protein_name} Ramachandran plot", subtitle, publication_stem)

    csv_path = data_dir / f"{safe_name}_phi_psi.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "protein": args.protein_name,
        "structure": str(args.structure.resolve()),
        "official_plot": str(args.official_plot.resolve()),
        "official_summary": str(args.summary.resolve()),
        "saves_job": stats["saves_job"],
        "total_residues_reported": stats["total_residues"],
        "phi_psi_pairs": len(rows),
        "non_gly_pro_denominator": denominator,
        "official_percentages": stats["percentages"],
        "official_counts_inferred_from_rounded_percentages": counts,
        "map_derived_counts": map_counts,
        "map_minus_official_count_discrepancy": discrepancies,
        "map_derived_disallowed_residues": disallowed_labels,
        "interpretation_boundary": "PROCHECK supports stereochemical plausibility but does not experimentally validate a predicted structure.",
    }
    (data_dir / f"{safe_name}_ramachandran_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_thesis_text(
        analysis_dir / f"毕业论文_{safe_name}_拉氏图_材料方法结果与图注.md",
        args.protein_name, stats, counts, denominator, disallowed_labels,
        args.saves_version, args.procheck_version,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
