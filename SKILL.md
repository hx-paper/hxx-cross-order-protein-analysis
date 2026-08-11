---
name: hxx-cross-order-protein-analysis
description: >-
  Analyze an insect protein from user-provided CDS and/or amino-acid sequence and produce HXX's standardized thesis-ready cross-order bioinformatics package, including NCBI homolog selection across Lepidoptera, Coleoptera, Diptera, Hymenoptera and Hemiptera, MUSCLE alignment, query-centred identity, IQ-TREE maximum-likelihood phylogeny with UFBoot, conserved-domain architecture, ProtParam-compatible physicochemical properties, compact ABC publication figures, three-line tables, and concise Chinese Methods, Results and legends. Use when the user says 分析我的蛋白, 蛋白跨目分析, 序列比对, 进化树加结构域, 理化性质, or asks to reproduce the finalized OBP/CREB/JHAMT/MRJP figure style.
---

# HXX cross-order protein analysis

Produce a complete, reproducible protein-analysis folder while keeping scientific claims conservative. Read [references/final-style-spec.md](references/final-style-spec.md) before drawing. Read [references/output-contract.md](references/output-contract.md) before creating files.

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

## 8. Draw and export

- Panel A: compact rectangular ML tree plus aligned domain architecture.
- Panel B: query-centred identity lollipop/dot plot grouped by insect order.
- Panel C: full-length alignment compressed to 2–3 blocks.
- Export separate A, B and C PDFs and one combined ABC figure in PDF, 300-dpi PNG and editable SVG. Also preserve the standalone tree-domain panel in all three formats.
- Apply the exact typography, palette, geometry and QA rules in [references/final-style-spec.md](references/final-style-spec.md).
- Use Matplotlib or equivalent vector-capable plotting. Embed editable text in PDF/SVG (`pdf.fonttype=42`, `svg.fonttype=none`) when possible.

## 9. Write thesis text

Write one concise Methods paragraph, one concise Results paragraph and one short figure legend. Include only the most informative values: protein length, identity range and maximum, key UFBoot value, domain coordinates and E-value. Avoid listing every order-specific range or weak deep node unless it changes the interpretation. Use “支持/表明/初步注释” rather than claiming verified biochemical function.

## 10. Verify before completion

- Render PDFs to images and inspect at final size.
- Confirm 10-pt body text is legible; bootstrap labels are 8 pt, transparent and off all lines; query is first, red and present in the domain panel; bars are aligned and not clipped; panels are compact; all names/accessions/domains match source data.
- Confirm all expected folders, raw data, scripts, separate PDFs, combined PDF/PNG/SVG, physicochemical TSV/JSON and concise thesis text exist.
- Report unresolved sequence, annotation or topology uncertainty explicitly.
