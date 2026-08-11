from __future__ import annotations

import re
import sys
from pathlib import Path


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: scaffold_protein_project.py <output-root> <protein-name>")
    root = Path(sys.argv[1]).resolve()
    protein = re.sub(r"[^A-Za-z0-9_.-]+", "_", sys.argv[2].strip())
    if not protein:
        raise SystemExit("Protein name is empty after sanitization")
    project = root / protein
    subdirs = [
        "analysis", "data/phylogeny", "data/smart_domains", "data/ibs2",
        "data/source_candidates", "data/structure_validation", "figures",
        "physicochemical_properties", "scripts",
    ]
    for subdir in subdirs:
        (project / subdir).mkdir(parents=True, exist_ok=True)
    print(project)


if __name__ == "__main__":
    main()
