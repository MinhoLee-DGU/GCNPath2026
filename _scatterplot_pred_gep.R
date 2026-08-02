#!/usr/bin/env Rscript

# ==============================================================================
# _scatterplot_pred_gep.R
# ==============================================================================
# GCNPath Utility: Scatter Plot of Gene/Pathway Expression vs. Predictions
#
# Plots gene or pathway expression values on the x-axis and model predictions
# on the y-axis, with support for filtering by subsets of cells or drugs,
# cell line annotation matching (TCGA Cancer Subtype), and target pathway coloring.
# ==============================================================================

# Suppress warnings and load libraries quietly
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
})

# Get script directory to resolve default annotation and mapping paths
initial_options = commandArgs(trailingOnly = FALSE)
file_arg_name = "--file="
script_match = grep(file_arg_name, initial_options)
if (length(script_match) > 0) {
  script_dir = dirname(sub(file_arg_name, "", initial_options[script_match]))
} else {
  script_dir = "."
}

# Helper function to print detailed usage instructions
print_usage = function() {
  cat("\n")
  cat("================================================================================\n")
  cat("                GCNPath Gene/Pathway GEX vs. Prediction Scatter Plot            \n")
  cat("================================================================================\n")
  cat("Usage:\n")
  cat("  Rscript _scatterplot_pred_gep.R --pred_file=<files> --gex_file=<file> --gex_name=<name> [options]\n\n")
  cat("Required Arguments:\n")
  cat("  --pred_file      Path to prediction file(s). Supports wildcards (e.g. *.csv)\n")
  cat("                   and multiple files separated by '|'.\n")
  cat("  --gex_file       Path to gene/pathway expression file (e.g. SANGER_RNA_TPM.csv).\n")
  cat("  --gex_name       Gene symbol (e.g. AKT1) or pathway name (e.g. BIOCARTA_AKT_PATHWAY).\n")
  cat("\n")
  cat("Optional Arguments:\n")
  cat("  --col_pred       Column name containing prediction values. (Default: 'Prediction')\n")
  cat("  --col_cell       Column name representing the cell/sample identifier. (Default: 'Cell')\n")
  cat("  --col_drug       Column name representing the drug identifier. (Default: 'Drug')\n")
  cat("  --subset_cell    Cells or Cancer Subtypes to select, separated by '|' or ','.\n")
  cat("  --subset_drug    Drugs or Target Pathways to select, separated by '|' or ','.\n")
  cat("  --anno_cells     Cell line annotation CSV path. (Default: NULL)\n")
  cat("  --col_anno_cells1 Cell ID column in cell annotation. (Default: 'SANGER_MODEL_ID')\n")
  cat("  --col_anno_cells2 Cancer type column in cell annotation. (Default: 'CANCER_TYPE')\n")
  cat("  --anno_drugs     Drug annotation CSV path. (Default: NULL)\n")
  cat("  --col_anno_drugs1 Drug ID column in drug annotation. (Default: 'Drug_CID')\n")
  cat("  --col_anno_drugs2 Target pathway column in drug annotation. (Default: 'Name')\n")
  cat("  --anno_genes     Gene annotation CSV path. (Default: 'data/cell/Anno_Genes.csv')\n")
  cat("  --col_anno_genes1 Gene symbol/ID column in gene annotation. (Default: 'ENTREZ_ID')\n")
  cat("  --col_anno_genes2 Entrez/Target ID column in gene annotation. (Default: 'HGNC_SYMBOL')\n")
  cat("  --out_file       Output plot path (e.g., .pdf, .png, .svg). (Default: 'scatterplot_pred_gep.pdf')\n")
  cat("  --width          Plot width in cm. (Default: 8)\n")
  cat("  --height         Plot height in cm. (Default: 8)\n")
  cat("  --title          Main title of the plot. (Default: NULL)\n")
  cat("  --xlab           X-axis label. (Default: gex_name)\n")
  cat("  --ylab           Y-axis label. (Default: 'Predicted ln(IC50)')\n")
  cat("  --fit_line       Add linear regression fit line (T/F). (Default: TRUE)\n")
  cat("  --show_stats     Show correlation statistics (PCC, SCC, p-value) (T/F). (Default: TRUE)\n")
  cat("  --stats_pos      Position of the stats box: 'topleft', 'topright', 'bottomleft', 'bottomright'. (Default: 'topleft')\n")
  cat("  --point_size     Size of scatter points. (Default: 1.5)\n")
  cat("  --point_alpha    Transparency of scatter points [0, 1]. (Default: 0.6)\n")
  cat("  --point_color    Point color or column name to color points dynamically. (Default: 'black')\n")
  cat("  --aggregate_mean Whether to aggregate predictions across seeds/folds by taking the mean\n")
  cat("                   grouped by cell and drug. (Default: TRUE)\n")
  cat("  --dpi            DPI for output raster plots. (Default: 500)\n")
  cat("  --rasterize      Rasterize points in SVG to reduce file size (T/F). (Default: TRUE)\n")
  cat("================================================================================\n\n")
}

# Parse command line arguments
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  print_usage()
  q(status = 0)
}

# Helper to parse arguments dynamically
parse_arg_value = function(arg_name, default_val, is_numeric = FALSE, is_logical = FALSE, allow_null = TRUE) {
  # Normalize dashes (replace Unicode en-dash, em-dash, and single hyphens with standard double hyphens for flags)
  clean_args = args
  clean_args = gsub("[\u2012\u2013\u2014\u2015]", "-", clean_args)
  clean_args = gsub("^-+", "--", clean_args)
  
  match_idx = which(startsWith(clean_args, paste0("--", arg_name, "=")))
  if (length(match_idx) > 0) {
    val = substring(args[match_idx[1]], regexpr("=", args[match_idx[1]]) + 1)
  } else {
    match_idx = which(clean_args == paste0("--", arg_name))
    if (length(match_idx) > 0 && match_idx[1] < length(args)) {
      val = args[match_idx[1] + 1]
    } else {
      return(default_val)
    }
  }
  
  if (allow_null && tolower(val) %in% c("null", "none", "")) return(NULL)
  if (is_numeric) return(as.numeric(val))
  if (is_logical) {
    if (tolower(val) %in% c("t", "true", "1", "yes")) return(TRUE)
    if (tolower(val) %in% c("f", "false", "0", "no")) return(FALSE)
    return(as.logical(val))
  }
  return(val)
}

# Retrieve argument values
pred_file       = parse_arg_value("pred_file", NULL)
gex_file        = parse_arg_value("gex_file", NULL)
gex_name        = parse_arg_value("gex_name", NULL)
col_pred        = parse_arg_value("col_pred", "Prediction")
col_cell        = parse_arg_value("col_cell", "Cell")
col_drug        = parse_arg_value("col_drug", "Drug")
subset_cell     = parse_arg_value("subset_cell", NULL)
subset_drug     = parse_arg_value("subset_drug", NULL)
anno_cells      = parse_arg_value("anno_cells", NULL)
col_anno_cells1 = parse_arg_value("col_anno_cells1", "SANGER_MODEL_ID")
col_anno_cells2 = parse_arg_value("col_anno_cells2", "CANCER_TYPE")
anno_drugs      = parse_arg_value("anno_drugs", NULL)
col_anno_drugs1 = parse_arg_value("col_anno_drugs1", "Drug_CID")
col_anno_drugs2 = parse_arg_value("col_anno_drugs2", "Name")
anno_genes      = parse_arg_value("anno_genes", "data/cell/Anno_Genes.csv")
col_anno_genes1 = parse_arg_value("col_anno_genes1", "ENTREZ_ID")
col_anno_genes2 = parse_arg_value("col_anno_genes2", "HGNC_SYMBOL")

out_file        = parse_arg_value("out_file", "scatterplot_pred_gep.pdf")
width           = parse_arg_value("width", 8, is_numeric = TRUE)
height          = parse_arg_value("height", 8, is_numeric = TRUE)
title           = parse_arg_value("title", NULL)
xlab            = parse_arg_value("xlab", NULL)
ylab            = parse_arg_value("ylab", "Predicted ln(IC50)")
fit_line        = parse_arg_value("fit_line", TRUE, is_logical = TRUE)
show_stats      = parse_arg_value("show_stats", TRUE, is_logical = TRUE)
stats_pos       = parse_arg_value("stats_pos", "topleft")
point_size      = parse_arg_value("point_size", 1.5, is_numeric = TRUE)
point_alpha     = parse_arg_value("point_alpha", 0.6, is_numeric = TRUE)
point_color     = parse_arg_value("point_color", "black")
aggregate_mean  = parse_arg_value("aggregate_mean", TRUE, is_logical = TRUE)
dpi             = parse_arg_value("dpi", 500, is_numeric = TRUE)
rasterize       = parse_arg_value("rasterize", TRUE, is_logical = TRUE)

# Validate required parameters
if (is.null(pred_file) || is.null(gex_file) || is.null(gex_name)) {
  cat("[Error] Missing required arguments: --pred_file, --gex_file, and --gex_name are mandatory.\n", file = stderr())
  print_usage()
  q(status = 1)
}

# Helper function to read file using fastest available package
read_csv_file = function(path, row_names = NULL) {
  if (!file.exists(path)) {
    stop(paste("[Error] File does not exist:", path))
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    df = data.table::fread(path, data.table = FALSE, header = TRUE)
  } else {
    df = read.csv(path, check.names = FALSE, header = TRUE)
  }
  
  if (!is.null(row_names) && row_names == 1) {
    rownames(df) = df[[1]]
    df = df[, -1, drop = FALSE]
  }
  return(df)
}

## 1. Load Expression Data
cat("Reading expression file:", gex_file, "...\n")
gex_df = read_csv_file(gex_file, row_names = 1)

# 2. Load Prediction Files (handles multiple paths, wildcards/globs, and regular expressions)
files_split = unlist(strsplit(pred_file, "\\|"))
files = c()
for (f in files_split) {
  # Literal path
  if (file.exists(f)) {
    files = c(files, f)
    next
  }
  
  # Sys.glob (standard wildcards like *, ?)
  globbed = Sys.glob(f)
  if (length(globbed) > 0) {
    files = c(files, globbed)
    next
  }
  
  # list.files with regular expression
  dir_name = dirname(f)
  pattern_name = basename(f)
  if (dir_name == "") dir_name = "."
  
  if (dir.exists(dir_name)) {
    # Match using user-provided regular expression directly
    matched = list.files(path = dir_name, pattern = pattern_name, full.names = TRUE)
    if (length(matched) > 0) {
      files = c(files, matched)
      next
    }
    
    # Try converting glob-style '*' to regex '.*'
    regex_pattern = gsub("\\*", ".*", pattern_name)
    matched_regex = list.files(path = dir_name, pattern = regex_pattern, full.names = TRUE)
    if (length(matched_regex) > 0) {
      files = c(files, matched_regex)
      next
    }
  }
  
  # Fallback to literal so we report missing file gracefully
  files = c(files, f)
}
files = unique(files)

# Exclude Grad-CAM or Integrated Gradient files from glob patterns if normal predictions are present
is_gcam_or_igrad = grepl("_gcam|_igrad|gcam|igrad", files)
if (any(!is_gcam_or_igrad)) {
  files = files[!is_gcam_or_igrad]
}

if (length(files) == 0) {
  stop(paste("[Error] No matching prediction files found for pattern/path:", pred_file))
}

cat("Found", length(files), "prediction file(s) to process.\n")

df_list = list()
for (idx in seq_along(files)) {
  f_path = files[idx]
  cat("Reading predictions from:", f_path, "...\n")
  temp_df = read_csv_file(f_path)
  
  # Validate cell, drug and prediction column existence
  if (!(col_cell %in% colnames(temp_df))) {
    stop(paste("[Error] Cell identifier column '", col_cell, "' not found in:", f_path))
  }
  if (!(col_drug %in% colnames(temp_df))) {
    stop(paste("[Error] Drug identifier column '", col_drug, "' not found in:", f_path))
  }
  if (!(col_pred %in% colnames(temp_df))) {
    stop(paste("[Error] Prediction column '", col_pred, "' not found in:", f_path))
  }
  
  df_list[[idx]] = temp_df
}

# Bind rows
pred_df = dplyr::bind_rows(df_list)

# Make sure identifiers are characters for robust matching
pred_df[[col_drug]] = as.character(pred_df[[col_drug]])
pred_df[[col_cell]] = as.character(pred_df[[col_cell]])

# Mean aggregate predictions over multiple seeds/folds if requested
if (aggregate_mean) {
  cat("Aggregating predictions across seeds/folds (taking mean)...\n")
  group_cols = c(col_cell, col_drug)
  
  pred_df = pred_df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      !!col_pred := mean(.data[[col_pred]], na.rm = TRUE),
      .groups = "drop"
    )
}

# 3. Load Annotations and Merge
# Load Cell Annotations
cell_anno_loaded = FALSE
if (!is.null(anno_cells) && file.exists(anno_cells)) {
  cat("Loading cell annotations from:", anno_cells, "...\n")
  cell_anno_df = read_csv_file(anno_cells)
  if (col_anno_cells1 %in% colnames(cell_anno_df)) {
    idx = match(pred_df[[col_cell]], cell_anno_df[[col_anno_cells1]])
    
    # Merge primary annotation column
    if (col_anno_cells2 %in% colnames(cell_anno_df)) {
      pred_df[[col_anno_cells2]] = cell_anno_df[[col_anno_cells2]][idx]
    }

    cell_anno_loaded = TRUE
  } else {
    warning(paste("Identifier column", col_anno_cells1, "not found in cell annotations."))
  }
} else {
  cat("Cell line annotations not found. Skipping merge.\n")
}

# Load Drug Annotations
drug_anno_loaded = FALSE
if (!is.null(anno_drugs) && file.exists(anno_drugs)) {
  cat("Loading drug annotations from:", anno_drugs, "...\n")
  drug_anno_df = read_csv_file(anno_drugs)
  # Ensure ID column is character for robust matching
  drug_anno_df[[col_anno_drugs1]] = as.character(drug_anno_df[[col_anno_drugs1]])
  
  if (col_anno_drugs1 %in% colnames(drug_anno_df)) {
    idx = match(pred_df[[col_drug]], drug_anno_df[[col_anno_drugs1]])
    
    # Merge primary annotation column
    if (col_anno_drugs2 %in% colnames(drug_anno_df)) {
      pred_df[[col_anno_drugs2]] = drug_anno_df[[col_anno_drugs2]][idx]
    }
    
    drug_anno_loaded = TRUE
  } else {
    warning(paste("Identifier column", col_anno_drugs1, "not found in drug annotations."))
  }
} else {
  cat("Drug annotations not found. Skipping merge.\n")
}

# 4. Filter by subset criteria (subset_cell / subset_drug)
if (!is.null(subset_cell) && subset_cell != "") {
  cells_selected = unlist(strsplit(subset_cell, "\\|"))
  cells_selected = unlist(strsplit(cells_selected, ","))
  cells_selected = trimws(cells_selected)
  cat("Subsetting cells using criteria:", paste(cells_selected, collapse = ", "), "\n")
  
  cell_mask = pred_df[[col_cell]] %in% cells_selected
  if (cell_anno_loaded) {
    if (col_anno_cells2 %in% colnames(pred_df)) {
      cell_mask = cell_mask | pred_df[[col_anno_cells2]] %in% cells_selected
    }
  }
  pred_df = pred_df[cell_mask, , drop = FALSE]
}

if (!is.null(subset_drug) && subset_drug != "") {
  drugs_selected = unlist(strsplit(subset_drug, "\\|"))
  drugs_selected = unlist(strsplit(drugs_selected, ","))
  drugs_selected = trimws(drugs_selected)
  cat("Subsetting drugs using criteria:", paste(drugs_selected, collapse = ", "), "\n")
  
  drug_mask = pred_df[[col_drug]] %in% drugs_selected
  if (drug_anno_loaded) {
    if (col_anno_drugs2 %in% colnames(pred_df)) {
      drug_mask = drug_mask | pred_df[[col_anno_drugs2]] %in% drugs_selected
    }
  }
  pred_df = pred_df[drug_mask, , drop = FALSE]
}

# Split gex_name by '|' to support multiple genes/identifiers
gex_names = unlist(strsplit(gex_name, "\\|"))
gex_names = trimws(gex_names)

cat("Processing GEX name(s):", paste(gex_names, collapse = ", "), "\n")

for (gex_name_single in gex_names) {
  cat("\n------------------------------------------------------------\n")
  cat("Processing:", gex_name_single, "\n")
  cat("------------------------------------------------------------\n")
  
  # Map GEX column (handling identifier mapping if columns do not match directly)
  mapped_col = NULL
  if (gex_name_single %in% colnames(gex_df)) {
    mapped_col = gex_name_single
  } else {
    # Try mapping from input name using anno_genes
    anno_genes_file = anno_genes
    
    # If the specified/default path does not exist, look up fallback paths
    if (is.null(anno_genes_file) || !file.exists(anno_genes_file)) {
      fallback_paths = c(
        "data/cell_data/Anno_Genes.csv",
        "../../processed_data/cell_data/SANGER_Passports/Anno_Genes.csv",
        "../processed_data/cell_data/SANGER_Passports/Anno_Genes.csv",
        "processed_data/cell_data/SANGER_Passports/Anno_Genes.csv"
      )
      for (p in fallback_paths) {
        if (file.exists(p)) {
          anno_genes_file = p
          break
        }
      }
    }
    
    if (!is.null(anno_genes_file) && file.exists(anno_genes_file)) {
      cat("Loading annotations from", anno_genes_file, "to map identifier:", gex_name_single, "...\n")
      anno_genes_df = read_csv_file(anno_genes_file)
      
      # Check if target columns exist
      c1 = col_anno_genes1
      c2 = col_anno_genes2
      
      # Try case-insensitive matching for column names if not exact match
      if (!(c1 %in% colnames(anno_genes_df))) {
        match_c1 = colnames(anno_genes_df)[tolower(colnames(anno_genes_df)) == tolower(c1)]
        if (length(match_c1) > 0) c1 = match_c1[1]
      }
      if (!(c2 %in% colnames(anno_genes_df))) {
        match_c2 = colnames(anno_genes_df)[tolower(colnames(anno_genes_df)) == tolower(c2)]
        if (length(match_c2) > 0) c2 = match_c2[1]
      }
      
      # If columns are found, perform mapping
      if (c1 %in% colnames(anno_genes_df) && c2 %in% colnames(anno_genes_df)) {
        match_row = which(tolower(as.character(anno_genes_df[[c2]])) == tolower(gex_name_single))
        if (length(match_row) > 0) {
          mapped_val = anno_genes_df[[c1]][match_row[1]]
          cat("Mapped", gex_name_single, "to:", mapped_val, "\n")
          if (as.character(mapped_val) %in% colnames(gex_df)) {
            mapped_col = as.character(mapped_val)
          }
        }
      } else {
        # Fallback to general lookup if column names are not matched
        match_col1 = NULL
        for (poss_col in c("HGNC_SYMBOL", "COSMIC_GENE_SYMBOL", "GENE_ID")) {
          if (poss_col %in% colnames(anno_genes_df)) {
            match_row = which(tolower(as.character(anno_genes_df[[poss_col]])) == tolower(gex_name_single))
            if (length(match_row) > 0) {
              match_col1 = poss_col
              break
            }
          }
        }
        
        match_col2 = if ("ENTREZ_ID" %in% colnames(anno_genes_df)) "ENTREZ_ID" else NULL
        if (is.null(match_col2)) {
          id_cols = colnames(anno_genes_df)[grepl("ENTREZ|ID", colnames(anno_genes_df), ignore.case = TRUE)]
          if (length(id_cols) > 0) match_col2 = id_cols[1]
        }
        
        if (!is.null(match_col1) && !is.null(match_col2)) {
          match_row = which(tolower(as.character(anno_genes_df[[match_col1]])) == tolower(gex_name_single))
          if (length(match_row) > 0) {
            mapped_val = anno_genes_df[[match_col2]][match_row[1]]
            cat("Mapped", gex_name_single, "using fallback columns (", match_col1, "->", match_col2, "):", mapped_val, "\n")
            if (as.character(mapped_val) %in% colnames(gex_df)) {
              mapped_col = as.character(mapped_val)
            }
          }
        }
      }
    }
  }
  
  # Case-insensitive column lookup if still not matched
  if (is.null(mapped_col)) {
    cols_lower = tolower(colnames(gex_df))
    match_idx = which(cols_lower == tolower(gex_name_single))
    if (length(match_idx) > 0) {
      mapped_col = colnames(gex_df)[match_idx[1]]
      cat("Found fallback column match:", mapped_col, "\n")
    }
  }
  
  if (is.null(mapped_col)) {
    warning(paste("[Warning] Column '", gex_name_single, "' not found in expression data. Skipping."))
    next
  }
  
  cat("Extracting GEX values for:", gex_name_single, "from column:", mapped_col, "\n")
  gex_vector = gex_df[[mapped_col]]
  names(gex_vector) = rownames(gex_df)

  # 5. Merge Expression Data
  pred_df_subset = pred_df # Copy to avoid altering baseline pred_df
  idx_gex = match(pred_df_subset[[col_cell]], names(gex_vector))
  pred_df_subset$Expression = gex_vector[idx_gex]
  
  # Clean up NAs
  plot_df = pred_df_subset %>%
    filter(!is.na(Expression) & !is.na(.data[[col_pred]]) & !is.infinite(Expression) & !is.infinite(.data[[col_pred]]))
  
  if (nrow(plot_df) == 0) {
    warning(paste("[Warning] No records left to plot for", gex_name_single, "after subsetting and matching. Skipping."))
    next
  }
  
  cat("Final dataset contains:", nrow(plot_df), "data points.\n")
  
  # 6. Calculate stats
  pcc_res = cor.test(plot_df$Expression, plot_df[[col_pred]], method = "pearson")
  scc_res = cor.test(plot_df$Expression, plot_df[[col_pred]], method = "spearman", exact = FALSE)
  
  pcc = pcc_res$estimate
  pcc_pval = pcc_res$p.value
  scc = scc_res$estimate
  
  cat(sprintf("Stats computed:\n  PCC = %.4f (p = %.3g)\n  SCC = %.4f\n", pcc, pcc_pval, scc))
  
  # 7. Determine Coloring Variable
  color_var = "none"
  if (!is.null(point_color) && point_color %in% colnames(plot_df)) {
    color_var = point_color
  }
  
  # Build plot mapping
  if (color_var != "none") {
    cat("Coloring points dynamically by:", color_var, "\n")
    plot_mapping = aes(x = Expression, y = .data[[col_pred]], color = .data[[color_var]])
  } else {
    cat("Coloring points with uniform color:", point_color, "\n")
    plot_mapping = aes(x = Expression, y = .data[[col_pred]])
  }
  
  # Determine current output file name for SVG check
  if (length(gex_names) > 1) {
    ext = tools::file_ext(out_file)
    base_name = tools::file_path_sans_ext(out_file)
    current_out_file_temp = paste0(base_name, "_", gex_name_single, ".", ext)
  } else {
    current_out_file_temp = out_file
  }
  
  save_svg = grepl("\\.svg$", current_out_file_temp, ignore.case = TRUE)
  use_raster = FALSE
  
  if (save_svg && nrow(plot_df) > 10000 && rasterize) {
    ggrastr_avail = tryCatch({
      requireNamespace("ggrastr", quietly = TRUE)
    }, error = function(e) {
      FALSE
    })
    
    if (ggrastr_avail) {
      use_raster = TRUE
      cat("SVG format requested with >10,000 points. Rasterizing the points layer using 'ggrastr' to optimize file size.\n")
    } else {
      warning("SVG format requested with >10,000 points, but the 'ggrastr' package is not fully loaded or installed. Points will be rendered as vector shapes.")
    }
  }

  pl = ggplot(plot_df, plot_mapping)
  
  # Add point layer (with optional rasterization)
  geom_point_layer = if (color_var != "none") {
    geom_point(alpha = point_alpha, size = point_size)
  } else {
    geom_point(alpha = point_alpha, size = point_size, color = point_color)
  }
  
  if (use_raster) {
    pl = pl + tryCatch({
      ggrastr::rasterise(geom_point_layer, dpi = dpi)
    }, error = function(e) {
      warning("Failed to rasterize with ggrastr due to loading error. Rendering standard vector points instead.")
      geom_point_layer
    })
  } else {
    pl = pl + geom_point_layer
  }
  
  # Add trend line
  if (fit_line) {
    pl = pl + geom_smooth(method = "lm", color = "#2980b9", fill = "#3498db", alpha = 0.15, linetype = "solid", linewidth = 1)
  }
  
  # Apply styling
  current_xlab = xlab
  if (is.null(current_xlab)) {
    current_xlab = gex_name_single
  }
  pl = pl + labs(
    title = title,
    x = current_xlab,
    y = ylab,
    color = if (color_var != "none") color_var else NULL
  )
  
  pl = pl + theme_pubr(legend = if (color_var != "none") "right" else "none") +
    theme(
      plot.title = element_text(face = "plain", size = 14, hjust = 0.5),
      axis.title = element_text(face = "plain", size = 12),
      axis.text = element_text(face = "plain", size = 10),
      legend.title = element_text(face = "plain", size = 10),
      legend.text = element_text(face = "plain", size = 9)
    )
  
  # Color scale if factor colored
  if (color_var != "none") {
    pl = pl + scale_color_brewer(palette = "Set2")
  }
  
  # Add statistics text box if requested
  if (show_stats) {
    min_x = min(plot_df$Expression)
    max_x = max(plot_df$Expression)
    min_y = min(plot_df[[col_pred]])
    max_y = max(plot_df[[col_pred]])
    
    if (stats_pos == "topleft") {
      x_pos = min_x + 0.05 * (max_x - min_x)
      y_pos = max_y - 0.05 * (max_y - min_y)
      hjust = 0
      vjust = 1
    } else if (stats_pos == "topright") {
      x_pos = max_x - 0.05 * (max_x - min_x)
      y_pos = max_y - 0.05 * (max_y - min_y)
      hjust = 1
      vjust = 1
    } else if (stats_pos == "bottomleft") {
      x_pos = min_x + 0.05 * (max_x - min_x)
      y_pos = min_y + 0.05 * (max_y - min_y)
      hjust = 0
      vjust = 0
    } else { # bottomright
      x_pos = max_x - 0.05 * (max_x - min_x)
      y_pos = min_y + 0.05 * (max_y - min_y)
      hjust = 1
      vjust = 0
    }
    
    annot_text = sprintf(
      "N = %s\nPCC = %.3f\nSCC = %.3f",
      format(nrow(plot_df), big.mark=","), pcc, scc
    )
    
    pl = pl + annotate(
      "label",
      x = x_pos, y = y_pos,
      label = annot_text,
      hjust = hjust, vjust = vjust,
      fill = "white", color = "black",
      alpha = 0.8,
      fontface = "plain", size = 3.5
    )
  }
  
  # Save output plot
  out_dir = dirname(out_file)
  if (out_dir != "." && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  if (length(gex_names) > 1) {
    ext = tools::file_ext(out_file)
    base_name = tools::file_path_sans_ext(out_file)
    current_out_file = paste0(base_name, "_", gex_name_single, ".", ext)
  } else {
    current_out_file = out_file
  }
  
  cat("Saving scatter plot to:", current_out_file, " [", width, "x", height, " cm ]\n")
  ggsave(
    filename = current_out_file,
    plot = pl,
    width = width,
    height = height,
    units = "cm",
    dpi = dpi
  )
}

cat("\nSuccess! Scatter plot(s) generated successfully.\n")
q(status = 0)
