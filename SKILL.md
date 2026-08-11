---
name: hxx-cross-order-protein-analysis
description: >-
  Analyze an insect protein from user-provided CDS, amino-acid sequence, and optional CIF/PDB structure and produce a standardized manuscript-ready bioinformatics package, including NCBI cross-order homolog selection, MUSCLE alignment, query-centred identity, IQ-TREE maximum-likelihood phylogeny with UFBoot, conserved-domain architecture, ProtParam-compatible physicochemical properties, SAVES/PROCHECK structure validation, classic red-yellow and publication teal-yellow Ramachandran plots, compact vector figures, three-line tables, and concise Chinese Methods, Results and legends for every figure. Use when the user says 分析我的蛋白, 蛋白跨目分析, 序列比对, 进化树加结构域, 理化性质, 拉氏图, Ramachandran plot, PROCHECK, 蛋白结构评价, or asks to reproduce the finalized Protein A/B/C/D anonymous figure style.
---

# Cross-order protein analysis

Produce a complete, reproducible protein-analysis folder while keeping scientific claims conservative. Read [references/final-style-spec.md](references/final-style-spec.md) before drawing, [references/ramachandran-workflow.md](references/ramachandran-workflow.md) before structure validation, [references/thesis-text-spec.md](references/thesis-text-spec.md) before writing, and [references/output-contract.md](references/output-contract.md) before creating files.

## 1. Validate input

- Accept CDS, amino-acid sequence, or both; remove whitespace and normalize case.
- When both are supplied, translate the CDS and compare it with the protein. Report mismatches, internal stops, frame problems, and the terminal stop codon before downstream work.
- Use the complete unmodified protein for alignment and physicochemical calculations unless the user explicitly requests mature-protein analysis.
- Preserve the user's target as `Query_<PROTEIN>_user_provided`; never silently replace it with a database sequence.

## 2. Build the project

Run `scripts/scaffold_protein_project.py <output-root> <protein-name>` and keep all later outputs inside that protein folder. Do not overwrite a confirmed final figure unless the user explicitly requests replacement; otherwise increment the version suffix.

## 3. Retrieve and select homologues

- Search NCBI Protein/RefSeq or BLASTP using the complete query sequence.
- Select 2–3 proteins from each of Lepidoptera, Coleoptera, Diptera, Hymenoptera and Hemiptera; use 3 per order when reliable homologues exist, normally 15 references plus the query.
- Prefer full-length, non-fragmented, well-annotated RefSeq/protein records with strong coverage and plausible length/domain architecture. Retain species, order, accession, description, length, coverage, E-value and identity in a TSV.
- Verify every accession and downloaded sequence against NCBI. Do not treat BLAST similarity alone as proof of orthology or identical function.
- For expanded families, screen for paralog mixing using annotations, domains, length and preliminary trees; describe the final tree as a family gene tree when appropriate.

## 4. Align and quantify similarity

- Align the query plus references with MUSCLE v5.3 using full-length amino-acid sequences.
- Calculate query-centred amino-acid identity only across positions where both sequences are non-gap; state this definition.
- Retain the full alignment for panel C. Compress long alignments into 2–3 horizontal blocks rather than displaying only a conserved fragment.
- Mark exact conservation when one non-gap residue occupies at least 70% of a column. Mark physicochemical similarity when one residue class occupies at least 70%. Use the established residue classes and visual rules from the style reference.
- Save the aligned FASTA, pairwise identity matrix, query-to-reference identity TSV and a JSON summary.

## 5. Infer phylogeny

- For tree inference, retain columns occupied in at least 70% of sequences unless the data justify another threshold; preserve and report the full and trimmed alignment sizes.
- Use IQ-TREE with ModelFinder/BIC and report the model actually selected. Run 1,000 ultrafast bootstrap replicates.
- Introduce only UFBoot in the thesis text and figure unless the user asks for SH-aLRT. Display only UFBoot values >=70.
- Draw a rectangular maximum-likelihood tree, not the obsolete circular version. Place the query as the first/top tip by rotating nodes without changing topology.
- Preserve branch-length meaning. Compact the display by shortening panel width and whitespace, not by falsifying distances. Include a substitutions/site scale bar.
- Put support values in 8 pt, with transparent background, offset away from branches and tip markers. Inspect every label manually, especially dense and bottom-most clades.

## 6. Annotate domains

- Submit every selected protein, including the query, to NCBI Batch CD-Search and consult CDD, SMART and Pfam evidence. Use SMART/Pfam family evidence where relevant, but do not print the database source beside every bar.
- Verify domain coordinates and names; use a clear full biological name in the legend rather than an opaque accession/code.
- Align domain tracks row-by-row with tree tips. Always show the query's domain track.
- Render wide glossy/gradient protein bars with restrained three-dimensional depth. Keep the domain panel close to the tree and about 25% shorter than early drafts while retaining proportional coordinates.

## 7. Calculate physicochemical properties

Run `scripts/calculate_physicochemical.py <protein.fasta> <protein-name> <output-dir>`. Report amino-acid length, molecular weight, theoretical pI, acidic/basic residues, molecular formula, atom count, 280-nm extinction coefficient, instability index, aliphatic index, GRAVY and N-end-rule half-life. Interpret instability index `<40` as predicted stable and `>40` as predicted unstable. State that calculations use the complete sequence without cleavage or post-translational modification.

## 8. Validate a predicted structure and draw Ramachandran plots

- When a CIF or PDB model is supplied, verify chain, residue count, missing backbone atoms and sequence correspondence before analysis. Do not infer a structure from sequence alone unless the user explicitly requests structure prediction.
- Convert AlphaFold/ModelCIF to standard PDB with `scripts/convert_modelcif_to_pdb.py <input.cif> <output.pdb>` before SAVES submission.
- Submit the PDB to UCLA SAVES and run PROCHECK. Preserve the original job URL, job number, version, official plot, summary and detailed report. Treat these files as the authoritative region statistics.
- Run `scripts/render_ramachandran_dual_palette.py` with the structure, official PROCHECK plot and summary. Generate both required renderings from the same official boundaries and statistics:
  1. classic PROCHECK-inspired red/yellow palette;
  2. publication teal/yellow palette.
- Export each rendering as 600-dpi PNG, PDF and editable SVG, plus residue-level phi/psi CSV, JSON summary and an automatically populated concise Chinese Methods/Results/legend file.
- Never present locally approximated favored/allowed regions as official PROCHECK scores. Explain that PROCHECK evaluates stereochemical plausibility and does not experimentally validate a predicted conformation. Follow [references/ramachandran-workflow.md](references/ramachandran-workflow.md).

## 9. Draw and export sequence figures

- Panel A: compact rectangular ML tree plus aligned domain architecture.
- Panel B: query-centred identity lollipop/dot plot grouped by insect order.
- Panel C: full-length alignment compressed to 2–3 blocks.
- Export separate A, B and C PDFs and one combined ABC figure in PDF, 300-dpi PNG and editable SVG. Also preserve the standalone tree-domain panel in all three formats.
- Apply the exact typography, palette, geometry and QA rules in [references/final-style-spec.md](references/final-style-spec.md).
- Use Matplotlib or equivalent vector-capable plotting. Embed editable text in PDF/SVG (`pdf.fonttype=42`, `svg.fonttype=none`) when possible.

## 10. Write thesis text for every output

Create one consolidated thesis file following [references/thesis-text-spec.md](references/thesis-text-spec.md). For every generated figure and physicochemical table, provide a short Methods paragraph, a short Results paragraph containing the actual values, and a figure/table legend. Cover tree-domain architecture, query-centred identity, full-length alignment, classic Ramachandran plot, publication Ramachandran plot and physicochemical table. Do not duplicate the same long explanation for the two Ramachandran palettes; state that they visualize identical data. Use “支持/表明/初步注释” rather than claiming verified biochemical function.

## 11. Verify before completion

- Render PDFs to images and inspect at final size.
- Confirm 10-pt body text is legible; bootstrap labels are 8 pt, transparent and off all lines; query is first, red and present in the domain panel; bars are aligned and not clipped; panels are compact; all names/accessions/domains match source data.
- Confirm all expected folders, raw data, scripts, separate PDFs, combined PDF/PNG/SVG, both Ramachandran palettes, official PROCHECK files, phi/psi CSV, physicochemical TSV/JSON and concise thesis text for every figure exist.
- Confirm the two Ramachandran versions have identical points, counts, flagged residues and percentages; only the palette and presentation may differ.
- Report unresolved sequence, annotation or topology uncertainty explicitly.
