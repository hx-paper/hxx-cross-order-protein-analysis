from __future__ import annotations

import argparse
import shlex
from pathlib import Path


def parse_atom_site(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    headers: list[str] = []
    start = None
    for i, line in enumerate(lines):
        if line.strip() == "_atom_site.group_PDB":
            j = i
            while j < len(lines) and lines[j].startswith("_atom_site."):
                headers.append(lines[j].strip().split(".", 1)[1])
                j += 1
            start = j
            break
    if start is None:
        raise ValueError("No _atom_site loop found in ModelCIF/mmCIF")
    atoms: list[dict[str, str]] = []
    for line in lines[start:]:
        stripped = line.strip()
        if not stripped or stripped == "#":
            if atoms:
                break
            continue
        if not (stripped.startswith("ATOM ") or stripped.startswith("HETATM ")):
            if atoms:
                break
            continue
        fields = shlex.split(stripped)
        if len(fields) != len(headers):
            raise ValueError(f"Unexpected atom_site row width: {len(fields)} != {len(headers)}")
        atoms.append(dict(zip(headers, fields)))
    if not atoms:
        raise ValueError("No coordinate records found")
    return atoms


def convert(source: Path, target: Path) -> tuple[int, int]:
    atoms = parse_atom_site(source)
    lines = [
        "REMARK   1 CONVERTED FROM MODELCIF/MMCIF FOR SAVES PROCHECK",
        "REMARK   2 B-FACTOR FIELD IS RETAINED FROM THE SOURCE MODEL",
    ]
    serial = 1
    last_chain = None
    last_resname = ""
    last_seq = 0
    residue_keys: set[tuple[str, int, str]] = set()
    for atom in atoms:
        if atom["group_PDB"] != "ATOM":
            continue
        chain = atom.get("auth_asym_id", atom.get("label_asym_id", "A"))[:1] or "A"
        resname = atom["label_comp_id"][:3]
        seq = int(atom.get("auth_seq_id", atom["label_seq_id"]))
        if last_chain is not None and chain != last_chain:
            lines.append(f"TER   {serial:5d}      {last_resname:>3s} {last_chain:1s}{last_seq:4d}")
            serial += 1
        name = atom["label_atom_id"]
        element = atom["type_symbol"].upper()[:2]
        atom_name = f" {name:<3s}" if len(name) < 4 and len(element) == 1 else f"{name:>4s}"
        altloc = " " if atom.get("label_alt_id", ".") in {".", "?"} else atom["label_alt_id"][:1]
        ins_value = atom.get("pdbx_PDB_ins_code", ".")
        icode = " " if ins_value in {".", "?"} else ins_value[:1]
        x, y, z = (float(atom[key]) for key in ("Cartn_x", "Cartn_y", "Cartn_z"))
        occupancy = float(atom.get("occupancy", 1.0))
        bfactor = float(atom.get("B_iso_or_equiv", 0.0))
        lines.append(
            f"ATOM  {serial:5d} {atom_name}{altloc}{resname:>3s} {chain:1s}{seq:4d}{icode}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          {element:>2s}  "
        )
        residue_keys.add((chain, seq, icode))
        serial += 1
        last_chain, last_resname, last_seq = chain, resname, seq
    if last_chain is not None:
        lines.append(f"TER   {serial:5d}      {last_resname:>3s} {last_chain:1s}{last_seq:4d}")
    lines.append("END")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines) + "\n", encoding="ascii")
    return serial - 1, len(residue_keys)


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert ModelCIF/mmCIF coordinates to a SAVES-compatible PDB")
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    args = parser.parse_args()
    records, residues = convert(args.source.resolve(), args.target.resolve())
    print(f"Wrote {args.target.resolve()} with {records} coordinate/TER records and {residues} residues")


if __name__ == "__main__":
    main()
