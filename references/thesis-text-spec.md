# Concise thesis text for all figures

Create `analysis/毕业论文_<PROTEIN>全部生物信息学图_材料方法结果与图注_简短版.md`. Use actual values from generated TSV/JSON/tree/domain/PROCHECK files. Do not leave invented placeholders in a delivered analysis.

## 1. Cross-order alignment and identity

- Methods: state NCBI source, five insect orders, 2–3 sequences per order, MUSCLE version, full-length alignment and pairwise identity definition.
- Results: give query length, identity range, highest species/accession/value and only one useful conservation statistic.
- Legend: identify identity and full-length-alignment panels and the conservation colors.

## 2. Phylogeny and domain architecture

- Methods: state occupancy trimming, IQ-TREE version, selected model, 1,000 UFBoot replicates, UFBoot display threshold, and CDD/SMART/Pfam domain annotation.
- Results: give the closest supported clade, its UFBoot, domain name/coordinates/E-value, and a one-sentence interpretation boundary.
- Legend: identify query highlighting, order colors, branch-length scale, UFBoot and domain tracks.

## 3. Physicochemical properties table

- Methods: state ExPASy ProtParam-compatible complete-sequence calculation without cleavage or post-translational modification.
- Results: report length, molecular weight, theoretical pI, instability index with stable/unstable interpretation, aliphatic index and GRAVY. Keep other values in the table.
- Table note: define the instability-index threshold and abbreviations.

## 4. Ramachandran structure validation

- Methods: state CIF-to-PDB conversion when applicable, SAVES/PROCHECK version and job, main-plot denominator, phi/psi recalculation, and dual-palette redraw using official boundaries.
- Results: report count and percentage in all four regions, name disallowed residues, and state whether outliers fall in low-confidence regions only when pLDDT/B-factor evidence exists.
- Classic legend: define red, yellow, pale yellow and white regions and residue marker shapes.
- Publication legend: define dark teal, light teal, pale yellow and white regions. State that it is an alternative visualization of the same official data, not a second analysis.
- Boundary: state that stereochemical quality does not replace experimental structure validation.

## Style

- Use one short paragraph for Methods and one short paragraph for Results per analysis block.
- Prefer 3–5 key numbers per paragraph. Put complete values in tables/legends.
- Introduce only UFBoot for the tree unless the user asks otherwise.
- Do not call a predicted function or conformation experimentally verified.
