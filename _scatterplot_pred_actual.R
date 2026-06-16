#!/usr/bin/env Rscript

# ==============================================================================
# plot_pred.R
# ==============================================================================
# GCNPath AI Model Prediction vs. True Scatter Plot Utility
#
# Designed for the Bio-protocol journal paper GCNPath protocol.
# Minimizes dependencies (uses only: ggplot2, ggpubr, dplyr).
# Supports flexible command-line arguments to handle any arbitrary input files,
# multi-fold cross-validation aggregation, dynamic color-coding, and custom dimensions.
# ==============================================================================

# Suppress warnings and load libraries quietly
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
})

# Helper function to print detailed usage instructions
print_usage = function() {
  cat("\n")
  cat("================================================================================\n")
  cat("                    GCNPath Prediction Scatter Plot Utility                     \n")
  cat("================================================================================\n")
  cat("This script generates a high-quality scatter plot comparing predicted and       \n")
  cat("actual drug sensitivity values (e.g., ln(IC50)) for the GCNPath AI model.       \n")
  cat("Designed specifically for publication-ready figures (e.g., Bio-protocol).\n\n")
  cat("Usage:\n")
  cat("  Rscript plot_pred.R --pred_file=<files> --col_pred=<col> --col_true=<col> [options]\n\n")
  cat("Required Arguments:\n")
  cat("  --pred_file      Path to prediction file(s). Supports wildcards (e.g. *.csv)   \n")
  cat("                   and multiple files separated by '|'.                          \n")
  cat("  --col_pred       Column name containing the predicted values.                  \n")
  cat("  --col_true       Column name containing the true/actual values.                \n")
  cat("\n")
  cat("Optional Arguments:\n")
  cat("  --col_cell       Column name representing the cell line identifier.            \n")
  cat("  --col_drug       Column name representing the drug identifier.                 \n")
  cat("  --aggregate_mean Whether to aggregate across folds/seeds by taking the mean    \n")
  cat("                   grouped by col_cell and col_drug. (Default: TRUE)             \n")
  cat("  --point_color    Point color. Either a single color name/hex code              \n")
  cat("                   (e.g., 'black', '#34495e') or a column name for dynamic      \n")
  cat("                   coloring. (Default: 'black')                                  \n")
  cat("  --point_size     Size of scatter points. (Default: 1.5)                        \n")
  cat("  --point_alpha    Transparency of scatter points [0, 1]. (Default: 0.25)        \n")
  cat("  --width          Plot width in cm. (Default: 10)                               \n")
  cat("  --height         Plot height in cm. (Default: 10)                              \n")
  cat("  --out_file       Output plot path. Supported: .pdf, .png, .svg, etc.           \n")
  cat("                   (Default: 'prediction_scatter_plot.pdf')                      \n")
  cat("  --title          Main title of the plot. (Default: NULL, no title plotted)     \n")
  cat("  --xlab           X-axis label. (Default: 'Actual ln(IC50)')                    \n")
  cat("  --ylab           Y-axis label. (Default: 'Predicted ln(IC50)')                 \n")
  cat("  --fit_line       Add linear regression fit line (T/F). (Default: FALSE)        \n")
  cat("  --xy_line        Add a diagonal reference line y = x (T/F). (Default: TRUE)    \n")
  cat("  --show_stats     Show performance metrics (PCC, SCC, RMSE) (T/F). (Default: TRUE)\n")
  cat("  --stats_pos      Position of the stats box. Options: 'topleft', 'topright',    \n")
  cat("                   'bottomleft', 'bottomright'. (Default: 'topleft')             \n")
  cat("================================================================================\n\n")
}

# Parse command line arguments
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  print_usage()
  q(status = 0)
}

# Helper to parse arguments dynamically
parse_arg_value = function(arg_name, default_val, is_numeric = FALSE, is_logical = FALSE) {
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
  
  if (is_numeric) return(as.numeric(val))
  if (is_logical) {
    if (tolower(val) %in% c("t", "true", "1", "yes")) return(TRUE)
    if (tolower(val) %in% c("f", "false", "0", "no")) return(FALSE)
    return(as.logical(val))
  }
  return(val)
}

# Set argument variables
pred_file = parse_arg_value("pred_file", NULL)
col_pred = parse_arg_value("col_pred", NULL)
col_true = parse_arg_value("col_true", NULL)
col_cell = parse_arg_value("col_cell", NULL)
col_drug = parse_arg_value("col_drug", NULL)
aggregate_mean = parse_arg_value("aggregate_mean", TRUE, is_logical = TRUE)
point_color = parse_arg_value("point_color", "black")
point_size = parse_arg_value("point_size", 1.5, is_numeric = TRUE)
point_alpha = parse_arg_value("point_alpha", 0.25, is_numeric = TRUE)
width = parse_arg_value("width", 10, is_numeric = TRUE)
height = parse_arg_value("height", 10, is_numeric = TRUE)
out_file = parse_arg_value("out_file", "prediction_scatter_plot.pdf")
dpi = parse_arg_value("dpi", 500, is_numeric = TRUE)
rasterize = parse_arg_value("rasterize", TRUE, is_logical = TRUE)
title = parse_arg_value("title", NULL)
xlab = parse_arg_value("xlab", "Actual ln(IC50)")
ylab = parse_arg_value("ylab", "Predicted ln(IC50)")
fit_line = parse_arg_value("fit_line", FALSE, is_logical = TRUE)
xy_line = parse_arg_value("xy_line", TRUE, is_logical = TRUE)
show_stats = parse_arg_value("show_stats", TRUE, is_logical = TRUE)
stats_pos = parse_arg_value("stats_pos", "topleft")

# Check required parameters
if (is.null(pred_file) || is.null(col_pred) || is.null(col_true)) {
  cat("[Error] Missing required arguments: --pred_file, --col_pred, and --col_true are mandatory.\n", file = stderr())
  print_usage()
  q(status = 1)
}

# 1. Load Prediction Files (handles multiple paths, wildcards/globs, and regular expressions)
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

if (length(files) == 0) {
  stop(paste("[Error] No matching files found for pattern/path:", pred_file))
}

cat("Found", length(files), "file(s) to process.\n")

# Reading utility using fastest available packages
read_file = function(path) {
  if (!file.exists(path)) {
    stop(paste("[Error] File does not exist:", path))
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, data.table = FALSE, header = TRUE))
  } else if (requireNamespace("readr", quietly = TRUE)) {
    return(as.data.frame(readr::read_csv(path, show_col_types = FALSE)))
  } else {
    return(read.csv(path, check.names = FALSE, header = TRUE))
  }
}

# Read and merge dataset
df_list = list()
for (idx in seq_along(files)) {
  f_path = files[idx]
  if (!file.exists(f_path)) {
    stop(paste("[Error] File does not exist:", f_path))
  }
  cat("Reading:", f_path, "...\n")
  temp_df = read_file(f_path)
  
  # Validate column existence
  if (!(col_pred %in% colnames(temp_df))) {
    stop(paste("[Error] Prediction column '", col_pred, "' not found in:", f_path))
  }
  if (!(col_true %in% colnames(temp_df))) {
    stop(paste("[Error] True value column '", col_true, "' not found in:", f_path))
  }
  
  # Auto-add fold/seed column if multiple files are loaded to prevent row conflicts
  if (length(files) > 1) {
    fold_col_name = "Fold_Temp_Auto_Index"
    if (!(fold_col_name %in% colnames(temp_df))) {
      # Extract number if any in basename, else use index
      seed_val = gsub("[^0-9]", "", basename(f_path))
      if (seed_val == "") {
        seed_val = as.character(idx)
      }
      temp_df[[fold_col_name]] = seed_val
    }
  }
  
  df_list[[idx]] = temp_df
}

# Bind rows safely using dplyr
df = dplyr::bind_rows(df_list)

# 2. Aggregation (Mean) over multiple seeds/folds grouped by Cell Line and Drug ID
if (aggregate_mean) {
  cell_present = !is.null(col_cell) && col_cell != "" && col_cell %in% colnames(df)
  drug_present = !is.null(col_drug) && col_drug != "" && col_drug %in% colnames(df)
  
  if (cell_present && drug_present) {
    cat(sprintf("Aggregating predictions across seeds/folds by taking mean grouped by Cell (%s) and Drug (%s)...\n", col_cell, col_drug))
    original_rows = nrow(df)
    
    # Identify grouping columns: col_cell, col_drug, and any other non-fold, non-seed metadata columns
    exclude_patterns = c("seed", "fold", "run", "index", "temp_auto")
    group_cols = colnames(df)
    group_cols = setdiff(group_cols, c(col_pred, col_true))
    
    # Filter out any dynamic fold/seed indices
    for (pat in exclude_patterns) {
      group_cols = group_cols[!grepl(pat, group_cols, ignore.case = TRUE)]
    }
    
    # Ensure always grouping by cell line and drug columns
    group_cols = unique(c(col_cell, col_drug, group_cols))
    group_cols = intersect(group_cols, colnames(df))
    
    # Perform aggregation
    df = df %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(
        !!col_pred := mean(.data[[col_pred]], na.rm = TRUE),
        !!col_true := mean(.data[[col_true]], na.rm = TRUE),
        .groups = "drop"
      )
    cat("Aggregation completed. Rows reduced from", original_rows, "to", nrow(df), ".\n")
  } else {
    if (length(files) > 1) {
      cat("[Info] Multiple files detected, but --col_cell and/or --col_drug was not specified. Plotting all predictions directly without aggregation.\n")
    }
  }
}

# 3. Clean Dataset of NAs and Infinite values
df_clean = df %>% 
  filter(!is.na(.data[[col_pred]]) & !is.infinite(.data[[col_pred]]) &
         !is.na(.data[[col_true]]) & !is.infinite(.data[[col_true]]))

if (nrow(df_clean) == 0) {
  stop("[Error] No valid rows left after filtering out NAs and infinite values.")
}

if (nrow(df) - nrow(df_clean) > 0) {
  cat("Filtered out", nrow(df) - nrow(df_clean), "rows containing NA/Inf values.\n")
}

# 4. Statistical Metrics Calculation
rmse = sqrt(mean((df_clean[[col_true]] - df_clean[[col_pred]])^2))
pcc = cor(df_clean[[col_true]], df_clean[[col_pred]], method = "pearson")
scc = cor(df_clean[[col_true]], df_clean[[col_pred]], method = "spearman")

cat(sprintf("\nMetrics calculated:\n  N = %s\n  RMSE = %.4f\n  PCC (Pearson) = %.4f\n  SCC (Spearman) = %.4f\n\n",
            format(nrow(df_clean), big.mark=","), rmse, pcc, scc))

# 5. Determine Point Color Mapping
is_color_column = FALSE
if (point_color %in% colnames(df_clean)) {
  is_color_column = TRUE
  cat("Coloring points dynamically using column:", point_color, "\n")
} else {
  cat("Coloring points with uniform color:", point_color, "\n")
}

# 6. Generate the Plot using ggplot2 & ggpubr
plot_mapping = if (is_color_column) {
  aes(x = .data[[col_true]], y = .data[[col_pred]], color = .data[[point_color]])
} else {
  aes(x = .data[[col_true]], y = .data[[col_pred]])
}

pl = ggplot(df_clean, plot_mapping)

# Check if output is SVG and should be rasterized for high points count to avoid massive vector files
save_svg = grepl("\\.svg$", out_file, ignore.case = TRUE)
use_raster = FALSE

if (save_svg && nrow(df_clean) > 10000 && rasterize) {
  ggrastr_avail = tryCatch({
    requireNamespace("ggrastr", quietly = TRUE)
  }, error = function(e) {
    FALSE
  })
  
  if (ggrastr_avail) {
    use_raster = TRUE
    cat("SVG format requested with >10,000 points. Rasterizing the points layer using 'ggrastr' to optimize file size.\n")
  } else {
    warning("SVG format requested with >10,000 points, but the 'ggrastr' package is not fully loaded or installed. Points will be rendered as vector shapes (output file size may be very large).")
  }
}

# Add point layer (with optional rasterization)
geom_point_layer = if (is_color_column) {
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

# Add diagonal y = x reference line (light red tomato color, thin, plain style)
if (xy_line) {
  pl = pl + geom_abline(slope = 1, intercept = 0, color = "red", alpha = 0.5, linewidth = 0.5, lty=2)
}

# Add linear regression fit line
if (fit_line) {
  pl = pl + geom_smooth(method = "lm", color = "blue", fill = "lightblue", se = TRUE, linetype = "solid", linewidth = 1)
}

# Apply custom styling & titles (plain, not bold)
pl = pl + labs(
  title = title,
  x = xlab,
  y = ylab
)

# Apply beautiful theme_pubr with 100% plain (not bold) font styling
pl = pl + theme_pubr(legend = ifelse(is_color_column, "right", "none")) +
  theme(
    plot.title = element_text(face = "plain", size = 14, hjust = 0.5),
    axis.title = element_text(face = "plain", size = 12),
    axis.text = element_text(face = "plain", size = 10)
  )

# Add professional color scales if coloring by column
if (is_color_column) {
  col_val = df_clean[[point_color]]
  if (is.character(col_val) || is.factor(col_val) || is.logical(col_val)) {
    pl = pl + scale_color_brewer(palette = "Set1")
  } else if (is.numeric(col_val)) {
    pl = pl + scale_color_viridis_c(option = "viridis")
  }
}

# Add statistics label box if requested (plain, not bold)
if (show_stats) {
  min_x = min(df_clean[[col_true]])
  max_x = max(df_clean[[col_true]])
  min_y = min(df_clean[[col_pred]])
  max_y = max(df_clean[[col_pred]])
  
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
    "N = %s\nPCC = %.3f\nSCC = %.3f\nRMSE = %.3f",
    format(nrow(df_clean), big.mark=","), pcc, scc, rmse
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

# 7. Save output plot
out_dir = dirname(out_file)
if (out_dir != "." && !dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("Saving plot to:", out_file, " [", width, "x", height, " cm ]\n")
ggsave(
  filename = out_file,
  plot = pl,
  width = width,
  height = height,
  units = "cm",
  dpi = dpi
)

cat("Success! Scatter plot generated successfully.\n")
