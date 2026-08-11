# Finalized figure style and revision ledger

## What changed during OBP, CREB, JHAMT and MRJP revision

- Increased all main labels and alignment text to 10 pt because the first combined figure was too small.
- Exported panels A, B and C separately as PDFs, plus a combined ABC figure.
- Replaced circular phylogenies in the combined figures with compact rectangular trees paired with domain architecture.
- Changed CREB from a short conserved-region display to a compressed full-length alignment, normally split across 2–3 blocks.
- Moved the target protein to the first/top position in every tree without changing topology.
- Reduced tree branch display width and outer whitespace repeatedly to make the whole figure compact.
- Moved domain tracks closer to the tree and shortened their display by about one quarter.
- Widened domain/protein bars and added restrained gradient/highlight/shadow effects for a polished IBS-like three-dimensional appearance.
- Replaced database codes in legends with simple, accurate biological domain names; retained full names when clarity required them.
- Removed repeated labels describing the domain-evaluation database from the figure itself.
- Restored the missing target-protein domain track and required it in every future figure.
- Reduced bootstrap text to 8 pt, removed white label boxes, and manually repositioned every number so none overlaps a branch. Dense CREB and bottom-clade regions require special inspection.
- Simplified thesis Methods and Results to short paragraphs with only key numerical evidence; mention UFBoot only.

## Fixed visual system

### Typography

- Font family: Arial for labels; DejaVu Sans Mono or another true monospaced font for alignment residues.
- Body, tip labels, axis labels and legends: 10 pt.
- Panel titles: 13 pt, bold, left aligned.
- Bootstrap labels: 8 pt, no box or background.
- Query label: 10 pt, bold, red.

### Insect-order palette

Use one palette in all panels and all proteins:

- Lepidoptera: `#E69F00`
- Coleoptera: `#0072B2`
- Diptera: `#009E73`
- Hymenoptera: `#CC79A7`
- Hemiptera: `#D55E00`
- Query: `#C62828`
- Branches: `#3F3F3F`

### Tree and domain panel

- Use a rectangular ML tree with thin dark-gray branches around 0.8–0.9 pt.
- Put the query at the top and optionally add a very pale red row highlight.
- Use small colored tip dots or order-color strips consistently.
- Italicize species names where the backend supports mixed typography; keep accessions upright.
- Show only UFBoot >=70 and place each number in clear space above/beside its node.
- Do not use a white rectangle behind support values.
- Keep the scale bar below the tree without excessive bottom whitespace.
- Align every structure bar precisely to its tree tip. Keep the tree-label-domain gaps narrow.
- Show protein length at the right edge in subtle gray when useful.
- Use a light gray glossy backbone and saturated glossy domain fill with top highlight and bottom shadow. Avoid cartoonish bevels.
- Keep bars visually wide enough to read and domain lengths proportional to the full protein.
- Use a short, accurate legend such as `MRJP/yellow family domain` or the full biological name `S-adenosyl-L-methionine-dependent methyltransferase domain`; never substitute a database accession for the name.

### Identity panel

- Use a clean horizontal lollipop/dot plot with light x-grid only.
- Group rows by insect order and color by the fixed palette.
- Abbreviate genus in crowded species labels; retain enough information to identify the accession in source data.
- Avoid decorative boxes unless an annotation is essential.

### Alignment panel

- Display the complete alignment. Split wide alignments into 2–3 compact blocks.
- Put the query first, bold and red-labeled.
- Use black background/white residue for exact conservation >=70%; medium gray background/dark residue for physicochemical similarity >=70%.
- Use orange for cysteine where biologically relevant; show signal peptide and mature-region cysteines for secreted OBPs when supported.
- Add subtle order-color strips and very light order-group shading; avoid heavy grid lines.

### Ramachandran plots

- Generate two standalone versions from identical phi/psi coordinates and official PROCHECK region statistics.
- Classic version: most favoured `#F50000`, additionally allowed `#FFF200`, generously allowed `#FFF9A6`, disallowed `#FFFFFF`.
- Publication version: most favoured `#5FA8A0`, additionally allowed `#B9DCD5`, generously allowed `#F4E6A2`, disallowed `#FAFAFA`.
- Use official PROCHECK/SAVES region boundaries. Do not substitute ellipses or locally approximated regions in thesis figures.
- Use Arial; axis and legend text about 10 pt, panel title 12–13 pt, and residue annotations 8 pt. Keep phi and psi axes square from -180° to 180°.
- Show other residues, glycine and proline with distinguishable shapes. Highlight and label disallowed residues; keep annotations off data points and axes.
- Add a horizontal four-class percentage bar below the plot. State the non-Gly/Pro denominator and the counts in most favoured/additionally allowed/generously allowed/disallowed order.
- Export each palette as 600-dpi PNG, PDF and editable SVG. Preserve the original official PROCHECK PNG/PDF separately.

### Deliverable QA

- Inspect PNG at 100% and render every PDF before delivery.
- Zoom into all bootstrap labels; no digit may touch or cross a vertical/horizontal branch.
- Inspect the bottom-most clade and scale bar separately.
- Check the target domain track, domain name, domain coordinates, protein lengths and row alignment.
- Confirm no obsolete circular tree remains in the final ABC figure.
- Compare both Ramachandran variants programmatically: point count, official percentages, class counts and flagged residues must match exactly.
