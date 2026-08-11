# Standard output contract

Create this structure under `<output-root>/<PROTEIN>/`:

```text
<PROTEIN>/
├── README.md
├── analysis/
│   ├── 毕业论文_<PROTEIN>组合图_材料方法结果与图注_简短版.md
│   └── 毕业论文_<PROTEIN>跨目序列及理化性质分析_材料方法与结果.docx
├── data/
│   ├── query_cds.fasta
│   ├── query_protein.fasta
│   ├── cross_order_<protein>_sequences.fasta
│   ├── cross_order_<protein>_alignment.fasta
│   ├── cross_order_<protein>_metadata.tsv
│   ├── query_to_reference_identity.tsv
│   ├── pairwise_identity_matrix.tsv
│   ├── selection_summary.json
│   ├── alignment_and_identity_summary.json
│   ├── phylogeny/
│   ├── smart_domains/
│   └── ibs2/
├── figures/
│   ├── Figure_A3_<PROTEIN>_ML_phylogeny_IBS2_domain_architecture_10pt.pdf
│   ├── Figure_A3_<PROTEIN>_ML_phylogeny_IBS2_domain_architecture_10pt.png
│   ├── Figure_A3_<PROTEIN>_ML_phylogeny_IBS2_domain_architecture_10pt.svg
│   ├── Figure_A_<PROTEIN>_ML_phylogeny_10pt.pdf
│   ├── Figure_B_<PROTEIN>_pairwise_identity_10pt.pdf
│   ├── Figure_C_<PROTEIN>_full_length_alignment_3line_10pt.pdf
│   ├── Figure_ABC_<PROTEIN>_cross_order_tree_domain_10pt.pdf
│   ├── Figure_ABC_<PROTEIN>_cross_order_tree_domain_10pt.png
│   └── Figure_ABC_<PROTEIN>_cross_order_tree_domain_10pt.svg
├── physicochemical_properties/
│   ├── <PROTEIN>_input_sequence.fasta
│   ├── <PROTEIN>_physicochemical_parameters.tsv
│   └── <PROTEIN>_physicochemical_results.json
└── scripts/
```

Use lowercase sanitized names inside data filenames and the user's canonical uppercase/locus form in visible labels. Add a version suffix instead of overwriting a final file unless replacement is explicitly authorized.
