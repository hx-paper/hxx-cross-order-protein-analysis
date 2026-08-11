# Ramachandran and PROCHECK workflow

## Evidence hierarchy

1. Treat the original SAVES/PROCHECK summary, detailed report and official plot as authoritative for region percentages and flags.
2. Recompute phi/psi angles from the supplied CIF/PDB to obtain editable scatter coordinates and a residue-level CSV.
3. Redraw the official categorical region background in two palettes. Keep points and official statistics identical between palettes.
4. Use local approximate favored/allowed regions only for exploratory checking. Never report their percentages as PROCHECK results.

## Required inputs

- Predicted or experimental structure in ModelCIF/mmCIF or PDB format.
- Official PROCHECK main Ramachandran PNG from SAVES.
- Official PROCHECK summary text. Preserve the SAVES job number and version when available.

If official SAVES files are not yet available, convert the structure to PDB and submit it. Do not fabricate official statistics while waiting.

## Conversion and submission

Run:

```text
python scripts/convert_modelcif_to_pdb.py model.cif model_for_SAVES.pdb
```

Check that chain identifiers, residue numbering, coordinate count and B-factor/model-confidence field are retained. Submit the PDB to UCLA SAVES, start PROCHECK, wait for completion, and download the main plot, summary and detailed report. Record the access date, SAVES version, PROCHECK version and job identifier.

## Dual-palette rendering

Run:

```text
python scripts/render_ramachandran_dual_palette.py STRUCTURE OFFICIAL_PNG SUMMARY_TXT FIGURE_DIR PROTEIN_NAME --data-dir STRUCTURE_VALIDATION_DIR --analysis-dir ANALYSIS_DIR
```

The renderer must:

- recover the four categorical boundaries from the official red/yellow plot;
- compute phi/psi without spanning chain breaks or missing peptide bonds;
- use the official summary percentages for the stacked bar;
- compare map-derived classifications with official rounded counts and report discrepancies;
- identify disallowed general residues from the official map and label them conservatively;
- export red/yellow and teal/yellow versions plus CSV, JSON and thesis text.

If the downloaded official PNG layout differs and the automatic plot crop fails, stop and inspect the crop rather than silently using an approximate background.

## Interpretation rules

- The main PROCHECK percentages use non-glycine and non-proline residues; state the denominator.
- Report all four classes: most favoured, additionally allowed, generously allowed and disallowed.
- Identify disallowed residues when the detailed report or map-based classification confirms them.
- The historical >90% most-favoured reference is a useful quality indicator, not a universal pass/fail rule for all predicted structures.
- For AlphaFold models, low pLDDT near an outlier may support a local-uncertainty interpretation, but does not prove why the geometry is unusual.
- A good Ramachandran distribution supports stereochemical plausibility only. It does not validate the biological state, oligomer, ligand pose or experimental accuracy.
