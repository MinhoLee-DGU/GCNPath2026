#!/usr/bin/env Rscript

# ==============================================================================
# barplot_importance.R
# ==============================================================================
# GCNPath Grad-CAM Pathway Importance Visualization Utility
#
# Generates high-quality, publication-ready horizontal barplots (and heatmaps)
# for Grad-CAM pathway importance scores. Supports dynamic scaling options 
# (raw, z-score, softmax, sum, min-max) and filtering by cell/drug.
# Outputs are saved in vector SVG format.
# ==============================================================================

# Suppress warnings and load libraries quietly
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(pheatmap)
})

# Helper function to print detailed usage instructions
print_usage = function() {
  cat("\n")
  cat("================================================================================\n")
  cat("                GCNPath Grad-CAM Importance Visualization Utility              \n")
  cat("================================================================================\n")
  cat("Generates horizontal barplots or heatmaps from Grad-CAM pathway importance CSVs.\n\n")
  cat("Usage:\n")
  cat("  Rscript barplot_importance.R --gcam_file=<files> [options]\n\n")
  cat("Required Arguments:\n")
  cat("  --gcam_file      Path to Grad-CAM file(s). Supports wildcards (e.g. *gcam.csv),\n")
  cat("                   and multiple files separated by '|'.\n")
  cat("\n")
  cat("Optional Arguments:\n")
  cat("  --col_cell       Column name representing the cell line ID. (Default: 'Cell')\n")
  cat("  --col_drug       Column name representing the drug ID. (Default: 'Drug')\n")
  cat("  --subset_cell    Cell line Name(s) or ID(s) to select, separated by '|'.\n")
  cat("                   (Default: NULL, no selection)\n")
  cat("  --subset_drug    Drug Name(s) or ID(s) to select, separated by '|'.\n")
  cat("                   (Default: NULL, no selection)\n")
  cat("  --anno_cells     Path to cell line annotation CSV for name mapping.\n")
  cat("                   (Default: '_case_study/Anno_Cells.csv')\n")
  cat("  --col_anno_cells1 Column name for Cell ID in cell annotation. (Default: 'SANGER_MODEL_ID')\n")
  cat("  --col_anno_cells2 Column name for Cell Name in cell annotation. (Default: 'MODEL_NAME')\n")
  cat("  --anno_drugs     Path to drug annotation CSV for name mapping.\n")
  cat("                   (Default: '_case_study/Anno_Drugs.csv')\n")
  cat("  --col_anno_drugs1 Column name for Drug ID in drug annotation. (Default: 'Drug_CID')\n")
  cat("  --col_anno_drugs2 Column name for Drug Name in drug annotation. (Default: 'Name')\n")
  cat("  --n_bar          Number of top pathways to select and display. Set to 0 to disable filtering. (Default: 10)\n")
  cat("  --scale          Scaling method for Grad-CAM scores. Since raw scores are very\n")
  cat("                   small, scaling improves readability. Options:\n")
  cat("                     'raw'     : No transformation (raw values).\n")
  cat("                     'zscore'  : standard z-score normalization.\n")
  cat("                     'softmax' : standard softmax transformation (sums to 1).\n")
  cat("                     'sum'     : simple relative proportion (sums to 1; recommended).\n")
  cat("                     'minmax'  : scale values to [0, 1] range.\n")
  cat("                   (Default: 'sum')\n")
  cat("  --heatmap        Generate a single heatmap (pheatmap) of union of top pathways\n")
  cat("                   instead of separate barplots when multiple cell-drug combinations\n")
  cat("                   are present (T/F). (Default: FALSE)\n")
  cat("  --combination    Check all Cartesian combinations of subset cells and drugs (T/F).\n")
  cat("                   If FALSE, 1-to-1 matching is used based on input order. (Default: FALSE)\n")
  cat("  --xlab_size      Font size of X-axis labels (Pathways) in heatmap. (Default: 9)\n")
  cat("  --ylab_size      Font size of Y-axis labels (Cell x Drug) in heatmap. (Default: 9)\n")
  cat("  --out_dir        Output directory for plots. (Default: 'plots_importance')\n")
  cat("  --width          Plot width in cm. (Default: 15)\n")
  cat("  --height         Plot height in cm. (Default: 10)\n")
  cat("  --dpi            DPI for raster plots (if converted). (Default: 500)\n")
  cat("  --svg            Save plots as SVG instead of PNG (T/F). (Default: FALSE)\n")
  cat("================================================================================\n\n")
}

# Get script directory to resolve relative paths
initial_options = commandArgs(trailingOnly = FALSE)
file_arg_name = "--file="
script_match = grep(file_arg_name, initial_options)
if (length(script_match) > 0) {
  script_dir = dirname(sub(file_arg_name, "", initial_options[script_match]))
} else {
  script_dir = "."
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

# Set argument variables
gcam_file = parse_arg_value("gcam_file", NULL, allow_null = FALSE)
col_cell = parse_arg_value("col_cell", "Cell")
col_drug = parse_arg_value("col_drug", "Drug")
subset_cell = parse_arg_value("subset_cell", NULL)
subset_drug = parse_arg_value("subset_drug", NULL)
anno_cells = parse_arg_value("anno_cells", "_case_study/Anno_Cells.csv")
col_anno_cells1 = parse_arg_value("col_anno_cells1", "SANGER_MODEL_ID")
col_anno_cells2 = parse_arg_value("col_anno_cells2", "MODEL_NAME")
anno_drugs = parse_arg_value("anno_drugs", "_case_study/Anno_Drugs.csv")
col_anno_drugs1 = parse_arg_value("col_anno_drugs1", "Drug_CID")
col_anno_drugs2 = parse_arg_value("col_anno_drugs2", "Name")
n_bar = parse_arg_value("n_bar", 10, is_numeric = TRUE)
scale_method = parse_arg_value("scale", "sum")
heatmap_opt = parse_arg_value("heatmap", FALSE, is_logical = TRUE)
combination = parse_arg_value("combination", FALSE, is_logical = TRUE)
xlab_size = parse_arg_value("xlab_size", 9, is_numeric = TRUE)
ylab_size = parse_arg_value("ylab_size", 9, is_numeric = TRUE)
out_dir = parse_arg_value("out_dir", "plots_importance")
width = parse_arg_value("width", 15, is_numeric = TRUE)
height = parse_arg_value("height", 10, is_numeric = TRUE)
dpi = parse_arg_value("dpi", 500, is_numeric = TRUE)
svg_opt = parse_arg_value("svg", FALSE, is_logical = TRUE)



# Ensure output directory exists
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# Helper to save plots
save_plot_safe = function(filename, plot, width, height, dpi) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "cm",
    dpi = dpi
  )
}

# Check required parameters
if (is.null(gcam_file)) {
  cat("[Error] Missing required argument: --gcam_file is mandatory.\n", file = stderr())
  print_usage()
  q(status = 1)
}

# 1. Load Grad-CAM Files (handles multiple paths, wildcards/globs, and regular expressions)
files_split = unlist(strsplit(gcam_file, "\\|"))
files = c()
for (f in files_split) {
  # Direct check
  if (file.exists(f)) {
    files = c(files, f)
    next
  }
  # If f is relative, check relative to script directory
  f_rel = file.path(script_dir, f)
  if (file.exists(f_rel)) {
    files = c(files, f_rel)
    next
  }
  
  # Try Sys.glob (standard wildcards like *, ?)
  globbed = Sys.glob(f)
  if (length(globbed) > 0) {
    files = c(files, globbed)
    next
  }
  globbed_rel = Sys.glob(f_rel)
  if (length(globbed_rel) > 0) {
    files = c(files, globbed_rel)
    next
  }
  
  # list.files with regular expression (direct and relative to script_dir)
  dir_name = dirname(f)
  pattern_name = basename(f)
  if (dir_name == "") dir_name = "."
  
  if (dir.exists(dir_name)) {
    matched = list.files(path = dir_name, pattern = pattern_name, full.names = TRUE)
    if (length(matched) > 0) {
      files = c(files, matched)
      next
    }
    regex_pattern = gsub("\\*", ".*", pattern_name)
    matched_regex = list.files(path = dir_name, pattern = regex_pattern, full.names = TRUE)
    if (length(matched_regex) > 0) {
      files = c(files, matched_regex)
      next
    }
  }
  
  dir_name_rel = dirname(f_rel)
  pattern_name_rel = basename(f_rel)
  if (dir_name_rel == "") dir_name_rel = "."
  
  if (dir.exists(dir_name_rel)) {
    matched_rel = list.files(path = dir_name_rel, pattern = pattern_name_rel, full.names = TRUE)
    if (length(matched_rel) > 0) {
      files = c(files, matched_rel)
      next
    }
    regex_pattern_rel = gsub("\\*", ".*", pattern_name_rel)
    matched_regex_rel = list.files(path = dir_name_rel, pattern = regex_pattern_rel, full.names = TRUE)
    if (length(matched_regex_rel) > 0) {
      files = c(files, matched_regex_rel)
      next
    }
  }
}

files = unique(files)
if (length(files) == 0) {
  cat(sprintf("[Error] No files found matching --gcam_file='%s'\n", gcam_file), file = stderr())
  q(status = 1)
}

cat(sprintf("[Info] Loading %d files...\n", length(files)))

# Load and aggregate across files/seeds
gcam_df = data.table()
for (i in 1:length(files)) {
  cat(sprintf("  [%d/%d] Loading: %s\n", i, length(files), basename(files[i])))
  temp = fread(files[i])
  gcam_df = rbind(gcam_df, temp)
}

# Aggregate by Cell and Drug (mean)
gcam_df = gcam_df[, lapply(.SD, mean), by = .(Cell, Drug)]

# 2. Annotations and Name Mapping
cat("[Info] Mapping Cell and Drug IDs using annotation files...\n")

# Initialize Cell_Name and Drug_Name with raw ID first
gcam_df$Cell_Name = as.character(gcam_df[[col_cell]])
gcam_df$Drug_Name = as.character(gcam_df[[col_drug]])

# Try cell mapping
if (!is.null(anno_cells) && file.exists(anno_cells)) {
  cell_ann = fread(anno_cells)
  if (col_anno_cells1 %in% colnames(cell_ann) && col_anno_cells2 %in% colnames(cell_ann)) {
    cell_map = setNames(as.character(cell_ann[[col_anno_cells2]]), as.character(cell_ann[[col_anno_cells1]]))
    mapped_cells = cell_map[gcam_df$Cell_Name]
    valid_mapped = !is.na(mapped_cells) & mapped_cells != "" & mapped_cells != "NA"
    gcam_df$Cell_Name[valid_mapped] = mapped_cells[valid_mapped]
  } else {
    cat(sprintf("[Warning] Columns '%s' or '%s' not found in %s. Using raw Cell IDs.\n", 
                col_anno_cells1, col_anno_cells2, anno_cells), file = stderr())
  }
} else {
  cat(sprintf("[Warning] Cell annotation file %s not found. Using raw Cell IDs.\n", anno_cells), file = stderr())
}

# Try file-based drug mapping
if (!is.null(anno_drugs) && file.exists(anno_drugs)) {
  drug_ann = fread(anno_drugs)
  if (col_anno_drugs1 %in% colnames(drug_ann) && col_anno_drugs2 %in% colnames(drug_ann)) {
    drug_map = setNames(as.character(drug_ann[[col_anno_drugs2]]), as.character(drug_ann[[col_anno_drugs1]]))
    file_mapped_drugs = drug_map[as.character(gcam_df[[col_drug]])]
    # File-based mapping should map any unmapped drug names
    valid_file_mapped = !is.na(file_mapped_drugs) & file_mapped_drugs != "" & file_mapped_drugs != "NA"
    gcam_df$Drug_Name[valid_file_mapped] = file_mapped_drugs[valid_file_mapped]
  } else {
    cat(sprintf("[Warning] Columns '%s' or '%s' not found in %s. Using raw Drug IDs.\n", 
                col_anno_drugs1, col_anno_drugs2, anno_drugs), file = stderr())
  }
} else {
  cat(sprintf("[Warning] Drug annotation file %s not found. Using raw Drug IDs.\n", anno_drugs), file = stderr())
}

# 3. Filtering by cell/drug subsets (supporting both '|' and ',' separators)
split_subset = function(val) {
  if (is.null(val)) return(NULL)
  unlist(strsplit(val, "[\\|,]+"))
}

cells_sel = split_subset(subset_cell)
drugs_sel = split_subset(subset_drug)

if (!combination && (!is.null(cells_sel) && !is.null(drugs_sel))) {
  # 1-to-1 matching based on input order (lengths may differ, truncate to shorter)
  len_cells = length(cells_sel)
  len_drugs = length(drugs_sel)
  if (len_cells != len_drugs) {
    min_len = min(len_cells, len_drugs)
    warning(sprintf("[Warning] subset_cell (length %d) and subset_drug (length %d) lengths do not match. Truncating to the shorter length: %d", len_cells, len_drugs, min_len))
    cells_sel = cells_sel[1:min_len]
    drugs_sel = drugs_sel[1:min_len]
  }
  
  keep_mask = rep(FALSE, nrow(gcam_df))
  for (i in seq_along(cells_sel)) {
    c_val = cells_sel[i]
    d_val = drugs_sel[i]
    
    cell_match = (gcam_df[[col_cell]] == c_val) | (gcam_df$Cell_Name == c_val)
    drug_match = (as.character(gcam_df[[col_drug]]) == d_val) | (gcam_df$Drug_Name == d_val)
    
    keep_mask = keep_mask | (cell_match & drug_match)
  }
  gcam_df = gcam_df[keep_mask, ]
} else {
  if (!combination) {
    warning("[Warning] combination=FALSE requires both --subset_cell and --subset_drug. Falling back to combination=TRUE behavior.")
  }
  if (!is.null(cells_sel)) {
    gcam_df = gcam_df %>% filter(.data[[col_cell]] %in% cells_sel | Cell_Name %in% cells_sel)
  }
  if (!is.null(drugs_sel)) {
    gcam_df = gcam_df %>% filter(as.character(.data[[col_drug]]) %in% drugs_sel | Drug_Name %in% drugs_sel)
  }
}

combinations = gcam_df %>% select(Cell = all_of(col_cell), Cell_Name, Drug = all_of(col_drug), Drug_Name) %>% distinct()

if (nrow(combinations) == 0) {
  cat("[Error] No matching Cell x Drug combinations found after filtering.\n", file = stderr())
  q(status = 1)
}

cat(sprintf("[Info] Found %d unique Cell x Drug combinations.\n", nrow(combinations)))

# Identify pathway columns (all columns except metadata)
metadata_cols = c(col_cell, col_drug, "Cell_Name", "Drug_Name", "Prediction", "LN_IC50", "Z_Pred", "Sample", "Subtype", "Status")
pathway_cols = setdiff(colnames(gcam_df), metadata_cols)

# Convert pathway columns to numeric matrix
gcam_mat = as.matrix(gcam_df[, pathway_cols, with = FALSE])

# 4. Apply Scaling per Combination
cat(sprintf("[Info] Applying '%s' scaling to pathway importance scores...\n", scale_method))
if (scale_method == "zscore") {
  scaled_mat = t(apply(gcam_mat, 1, function(x) {
    s = sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
  }))
} else if (scale_method == "softmax") {
  scaled_mat = t(apply(gcam_mat, 1, function(x) {
    ex = exp(x - max(x, na.rm = TRUE))
    ex / sum(ex, na.rm = TRUE)
  }))
} else if (scale_method == "sum") {
  scaled_mat = t(apply(gcam_mat, 1, function(x) {
    s = sum(x, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(x)) else x / s
  }))
} else if (scale_method == "minmax") {
  scaled_mat = t(apply(gcam_mat, 1, function(x) {
    r = diff(range(x, na.rm = TRUE))
    if (is.na(r) || r == 0) rep(0, length(x)) else (x - min(x, na.rm = TRUE)) / r
  }))
} else { # raw
  scaled_mat = gcam_mat
}

scaled_gcam_df = copy(gcam_df)
scaled_gcam_df[, (pathway_cols) := as.data.table(scaled_mat)]

# 5. Plotting Logic
if (!heatmap_opt) {
  # Scenario A: Draw separate barplots for each combination (Output: SVG, Plain Fonts, Added Margins)
  for (i in 1:nrow(combinations)) {
    comb = combinations[i, ]
    cat(sprintf("[Info] Plotting Barplot for: %s x %s\n", comb$Cell_Name, comb$Drug_Name))
    
    row_data = scaled_gcam_df %>% filter(.data[[col_cell]] == comb$Cell & .data[[col_drug]] == comb$Drug)
    values = unlist(row_data[, pathway_cols, with = FALSE])
    names(values) = pathway_cols
    
    # Filter for values > 0 (keep only positive ones)
    valid_idx = which(values > 0)
    if (length(valid_idx) == 0) {
      cat(sprintf("  [Warning] No positive pathways for %s x %s. Skipping.\n", comb$Cell_Name, comb$Drug_Name))
      next
    }
    
    values = values[valid_idx]
    sorted_values = sort(values, decreasing = TRUE)
    top_values = if (n_bar <= 0) sorted_values else head(sorted_values, n_bar)
    
    plot_df = data.frame(
      Pathway = names(top_values),
      Value = as.numeric(top_values),
      stringsAsFactors = FALSE
    )
    
    # Ensure descending order from top to bottom
    plot_df$Pathway = factor(plot_df$Pathway, levels = rev(plot_df$Pathway))
    
    p = ggplot(plot_df, aes(x = Pathway, y = Value)) +
      geom_bar(stat = "identity", fill = "#3b82f6", width = 0.7, alpha = 0.9) +
      coord_flip() +
      theme_minimal(base_size = 11) +
      theme(
        # Set all text elements to plain font style globally to guarantee plain font
        text = element_text(face = "plain"),
        axis.text.y = element_text(color = "grey10", size = 9, face = "plain", margin = margin(r = 5)),
        axis.text.x = element_text(color = "grey10", size = 9, face = "plain"),
        # Generous top margins for both x and y titles (for original and flipped axes)
        # to push the "Grad-CAM Score" label down and create elegant breathing room/margin
        axis.title.x = element_text(margin = margin(t = 15), face = "plain"),
        axis.title.y = element_text(margin = margin(t = 15, r = 15), face = "plain"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain", margin = margin(b = 10)),
        plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5, face = "plain", margin = margin(b = 10))
      ) +
      labs(
        x = NULL, 
        y = paste0("Grad-CAM Score (Scale: ", scale_method, ")"),
        title = paste0("Grad-CAM Pathway Importance"),
        subtitle = paste0("Cell: ", comb$Cell_Name, " | Drug: ", comb$Drug_Name)
      )
    
    sanitized_cell = gsub("[^A-Za-z0-9_-]", "_", comb$Cell_Name)
    sanitized_drug = gsub("[^A-Za-z0-9_-]", "_", comb$Drug_Name)
    ext = if (svg_opt) "svg" else "png"
    filename = file.path(out_dir, sprintf("barplot_%s_%s.%s", sanitized_cell, sanitized_drug, ext))
    
    save_plot_safe(filename, p, width = width, height = height, dpi = dpi)
    cat(sprintf("  [Success] Saved to %s\n", filename))
  }
} else {
  # Scenario B: Draw a single pheatmap of the union of top pathways (Output: SVG, Transposed with Cell x Drug in Rows)
  cat("[Info] Plotting Heatmap with Cell x Drug in Rows...\n")
  all_top_paths = c()
  comb_names = c()
  comb_list = list()
  
  for (i in 1:nrow(combinations)) {
    comb = combinations[i, ]
    comb_name = paste0(comb$Cell_Name, " x ", comb$Drug_Name)
    comb_names = c(comb_names, comb_name)
    
    row_data = scaled_gcam_df %>% filter(.data[[col_cell]] == comb$Cell & .data[[col_drug]] == comb$Drug)
    values = unlist(row_data[, pathway_cols, with = FALSE])
    names(values) = pathway_cols
    
    valid_values = values[values > 0]
    sorted_values = sort(valid_values, decreasing = TRUE)
    top_paths = if (n_bar <= 0) names(sorted_values) else names(head(sorted_values, n_bar))
    
    all_top_paths = c(all_top_paths, top_paths)
    comb_list[[comb_name]] = values
  }
  
  union_paths = unique(all_top_paths)
  
  if (length(union_paths) == 0) {
    cat("[Error] No positive pathways found across selected combinations. Heatmap aborted.\n", file = stderr())
    q(status = 1)
  }
  
  cat(sprintf("[Info] Union of top %d pathways covers %d unique pathways.\n", n_bar, length(union_paths)))
  
  # Build heatmap matrix (Rows: Cell x Drug combinations, Columns: Pathways)
  heat_mat = matrix(0, nrow = length(comb_names), ncol = length(union_paths))
  rownames(heat_mat) = comb_names
  colnames(heat_mat) = union_paths
  
  for (comb_name in comb_names) {
    heat_mat[comb_name, ] = comb_list[[comb_name]][union_paths]
  }
  
  heat_mat[is.na(heat_mat)] = 0
  
  ext = if (svg_opt) "svg" else "png"
  filename = file.path(out_dir, paste0("heatmap_importance.", ext))
  
  # Set up the appropriate graphics device
  if (svg_opt) {
    svg(filename, width = width / 2.54, height = height / 2.54)
  } else {
    png(filename, width = width / 2.54, height = height / 2.54, units = "in", res = dpi)
  }
  
  # Use main with trailing newlines (\n\n) to create vertical spacing/margin below title
  title_text = paste0("Pathway Importance Heatmap (Scale: ", scale_method, ")\n\n")
  
  pheatmap(
    heat_mat,
    main = title_text,
    color = colorRampPalette(c("white", "palegoldenrod", "firebrick3"))(100),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    angle_col = 45,
    display_numbers = FALSE,
    fontsize = 9,
    fontsize_row = ylab_size, # Rows correspond to Y-axis (Cell x Drug)
    fontsize_col = xlab_size  # Columns correspond to X-axis (Pathways)
  )
  dev.off()
  
  cat(sprintf("[Success] Saved heatmap to %s\n", filename))
}

cat("[Info] All operations completed successfully!\n")
