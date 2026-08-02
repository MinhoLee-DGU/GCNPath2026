# GCNPath Gene/Pathway Expression vs. Prediction Scatter Plot Manual (`_scatterplot_pred_gep.R`)

Relate predicted response values (e.g. ln(IC50)) to molecular features such as gene expression level (TPM) or pathway activation scores (GSVA). This utility maps identifiers using gene, cell, and drug annotation matrices, filters combinations by cell line or drug subsets, aggregates predictions, and generates publication-ready scatter plots.

---

## 💻 CLI Usage Examples

### Example A: Relate Predicted Cisplatin Sensitivity to AKT1 Gene Expression (Figure 6A)
Evaluate how SCLC cell lines respond to cisplatin compared to their level of `AKT1` gene expression:

```bash
Rscript _scatterplot_pred_gep.R \
  --pred_file="results/IC50_GDSC/Normal/RGCN/pred_sclc_seed[0-9]+.csv" \
  --gex_file="data/cell_data/SANGER_RNA_TPM_Filt.csv" \
  --gex_name="AKT1" \
  --subset_drug="Cisplatin" \
  --anno_drugs="data/drug_data/Anno_Drugs.csv" \
  --col_anno_drugs1="Drug_CID" \
  --col_anno_drugs2="Name" \
  --anno_genes="data/cell_data/Anno_Genes.csv" \
  --col_anno_genes1="ENTREZ_ID" \
  --col_anno_genes2="HGNC_SYMBOL" \
  --out_file="pred_cisplatin_akt1.png"
```

*Note:*
- In the prediction CSV, drugs are identified as PubChem CIDs (e.g. `5702198`). Using `--anno_drugs` mapping, you can specify `"Cisplatin"` directly in `--subset_drug`.
- In the gene expression file, genes are identified as Entrez IDs (e.g. `207`). Using `--anno_genes` mapping, you can specify `"AKT1"` directly in `--gex_name`.

### Example B: Relate Predicted Cisplatin Sensitivity to BIOCARTA AKT Pathway Score (Figure 6B)
Evaluate predicted cisplatin sensitivity against BIOCARTA AKT pathway activation scores:

```bash
Rscript _scatterplot_pred_gep.R \
  --pred_file="results/IC50_GDSC/Normal/RGCN/pred_sclc_seed[0-9]+.csv" \
  --gex_file="processed/cell_data_biocarta/SANGER_RNA_GSVA.csv" \
  --gex_name="BIOCARTA_AKT_PATHWAY" \
  --subset_drug="Cisplatin" \
  --anno_drugs="data/drug_data/Anno_Drugs.csv" \
  --col_anno_drugs1="Drug_CID" \
  --col_anno_drugs2="Name" \
  --out_file="pred_cisplatin_akt_pathway.png"
```

*Note:*
- Since the gene expression/pathway score file (`SANGER_RNA_GSVA.csv`) uses standard names matching `--gex_name` directly, the `--anno_genes` arguments are omitted.

---

## 📥 Input & 📤 Output Description

```text
Rscript _scatterplot_pred_gep.R \
  --pred_file      [Path to Prediction CSV file(s); Supports wildcards and '|' separators; Mandatory] \
  --gex_file       [Path to Gene/Pathway Expression CSV file; Mandatory] \
  --gex_name       [Target Gene Symbol or Pathway Name to plot on X-Axis; Mandatory] \
  --col_pred       [Column Name of predicted values in --pred_file; Default: 'Prediction'] \
  --col_cell       [Column Name of cell line IDs in --pred_file; Default: 'Cell'] \
  --col_drug       [Column Name of drug IDs in --pred_file; Default: 'Drug'] \
  --subset_cell    [Cell Line ID(s) or Cancer Subtypes to select, separated by '|'] \
  --subset_drug    [Drug ID(s) or Target Pathways to select, separated by '|'] \
  --anno_cells     [Path to Cell Line Annotation CSV file; Default: NULL] \
  --col_anno_cells1 [Cell ID column name in cell annotations; Default: 'SANGER_MODEL_ID'] \
  --col_anno_cells2 [Cancer type/name column name in cell annotations; Default: 'CANCER_TYPE'] \
  --anno_drugs     [Path to Drug Annotation CSV file; Default: NULL] \
  --col_anno_drugs1 [Drug ID column name in drug annotations; Default: 'Drug_CID'] \
  --col_anno_drugs2 [Drug Name column name in drug annotations; Default: 'Name'] \
  --anno_genes     [Path to Gene Annotation CSV file; Default: 'data/cell/Anno_Genes.csv'] \
  --col_anno_genes1 [Gene Symbol/ID column in gene annotations; Default: 'ENTREZ_ID'] \
  --col_anno_genes2 [Entrez/Target ID column in gene annotations; Default: 'HGNC_SYMBOL'] \
  --out_file       [Output Plot Path; Default: 'scatterplot_pred_gep.pdf'] \
  --width          [Plot Width in cm; Default: 8] \
  --height         [Plot Height in cm; Default: 8] \
  --title          [Main Title of the Plot; Default: NULL] \
  --xlab           [X-axis Label; Default: value of --gex_name] \
  --ylab           [Y-axis Label; Default: 'Predicted ln(IC50)'] \
  --fit_line       [Add Linear Regression Fit Line (T/F); Default: TRUE] \
  --show_stats     [Show Correlation Stats (PCC, SCC, p-value) (T/F); Default: TRUE] \
  --stats_pos      [Position of Stats Box ('topleft', 'topright', 'bottomleft', 'bottomright'); Default: 'topleft'] \
  --point_size     [Size of Scatter Points; Default: 1.5] \
  --point_alpha    [Transparency of Scatter Points [0, 1]; Default: 0.6] \
  --point_color    [Point Color or column name to color points dynamically; Default: 'black'] \
  --aggregate_mean [Whether to Aggregate predictions across seeds/folds (T/F); Default: TRUE] \
  --dpi            [Plot DPI; Default: 500] \
  --rasterize      [Rasterize Points layer in SVG outputs (T/F); Default: TRUE]
```
