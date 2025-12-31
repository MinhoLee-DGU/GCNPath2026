#!/usr/bin/env Rscript

##### 1. Preparation

suppressMessages(library(igraph))
suppressMessages(library(ggpubr))
suppressMessages(library(cogena))
suppressMessages(library(graphite))
suppressMessages(library(reshape2))

source("../functions.R")
source("functions.R")
loadings()


dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/TPM_Ent.csv", dir)
SANGER_RNA = fread_def(file, check_names=F, header=T)

dir = "../../processed_data/cell_data/CCLE_DepMap"
file = sprintf("%s/TPM_Ent.csv", dir)
CCLE_RNA = fread_def(file, check_names=F, header=T)

dir = "../../processed_data/cell_data/GDSC"
file = sprintf("%s/RNA_Array_Ent.csv", dir)
GDSC_RNA = fread_def(file, check_names=F, header=T)

dir = "../../raw_data/MSigDB"
file = sprintf("%s/c2.cp.biocarta.v2023.1.Hs.entrez.gmt", dir)
Path_List = gmt2list(file)
path_genes = Path_List %>% unlist %>% unique   # 1509

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)

# dir = "../../processed_data/cell_data/SANGER_Passports"
# file = sprintf("%s/Anno_Cells.csv", dir)
# Anno_Cells = read.csv(file)


# Examine ArcSinh normalization adopted in PaccMann and PaccMann_SG
# See project/_prepare_paccmann_sanger/code.R
# https://github.com/PaccMann/paccmann_predictor/issues/6

examine_asinh = T
if (examine_asinh) {
  SANGER_RNA_Asinh = 2**SANGER_RNA-1
  SANGER_RNA_Asinh %>% rowSums %>% range %>% round(3)   # 969858.9 1000014.3
  SANGER_RNA_Asinh = SANGER_RNA_Asinh %>% asinh %>% as.data.frame
  SANGER_RNA_Asinh %>% rowSums %>% range %>% round(3)   # 25700.88 62454.64
  
  CCLE_RNA_Asinh = 2**CCLE_RNA-1
  CCLE_RNA_Asinh %>% rowSums %>% range %>% round(3)   # 652927.2 960343.3
  CCLE_RNA_Asinh = CCLE_RNA_Asinh %>% asinh %>% as.data.frame
  CCLE_RNA_Asinh %>% rowSums %>% range %>% round(3)   # 17002.57 50365.57
}


### Define statistic functions
# Those functions are already defined in 0_functions.R

# RMSE_Norm [NRMSE-IQR]
# RNSE Normalized by the interquartile range [IQR] of observations

# stat_pair
# Apply a statistic function (such as RMSE, R2, Euclidean distance)
# to all rows or columns of two data.frames in pairwise manner

# stat_pair_apply
# Apply a statistic function (such as RMSE, R2, Euclidean distance)
# to all rows or columns of two data.frames in matched rows or columns

# stat_self
# Apply a statistic function (such as RMSE, R2, Euclidean distance)
# to all rows or columns of one data.frame in pairwise manner

# corr_self
# Apply a correlation function (method="pearson", use="pairwise.complete.obs")
# to all rows or columns of one data.frame in pairwise manner



##### 2. Calculate pathway activity scores [GSVA]

### GSVA vs ssGSEA & singscore
rna_genes = SANGER_RNA %>% colnames
sum(path_genes %in% rna_genes)   # 99.80% [1506/1509]

Path_List %>% sapply(length) %>% range
Path_List %>% sapply(function(x) length(intersect(x, rna_genes))) %>% range
Path_List %>% sapply(function(x) length(intersect(x, rna_genes))/length(x)) %>% hist
# [5, 81], [5, 81], Mostly above 0.95


### 2-1. GSVA is sensitive to the number of total genes in omics data? [X]

# Define stable genes as those with lowest sd
stable = getStableGenes(n_stable=1000, type="carcinoma")   # 1000
idx = match(stable, Anno_Genes$HGNC_SYMBOL)
stable = Anno_Genes$ENTREZ_ID[idx] %>% na.omit   # 1000

sum(stable %in% colnames(SANGER_RNA))   # 1000
sum(stable %in% colnames(GDSC_RNA))     # 953
sum(stable %in% colnames(CCLE_RNA))     # 999
sum(stable %in% path_genes)             # 119

stable_100 = stable %>% 
  intersect(path_genes) %>% 
  intersect(colnames(SANGER_RNA)) %>% head(100)   # 100


# cf. ssGSEA scores are reproducible, regardless of sample composition? [X. O if ssgsea.norm=F]
# ssGSEA = SANGER_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea")
# ssGSEA_100 = SANGER_RNA[1:100, ] %>% gsva_def(Path_List, filt_genes=T, method="ssgsea")
# identical(ssGSEA[1:100, ], ssGSEA_100)   # F
# 
# ssGSEA[1:3, 1:2]
# BIOCARTA_GRANULOCYTES_PATHWAY BIOCARTA_LYM_PATHWAY 
# SIDM00001                    0.07075104           0.08435927       
# SIDM00002                   -0.25323122          -0.15020834
# SIDM00003                   -0.27420446          -0.16196303
# 
# ssGSEA_100[1:3, 1:2]
# BIOCARTA_GRANULOCYTES_PATHWAY BIOCARTA_LYM_PATHWAY 
# SIDM00001                    0.07535252            0.0898458
# SIDM00002                   -0.26970079           -0.1599775
# SIDM00003                   -0.29203808           -0.1724967
#
# ssGSEA = SANGER_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea", ssgsea.norm=F)
# ssGSEA_100 = SANGER_RNA[1:100, ] %>% gsva_def(Path_List, filt_genes=T, method="ssgsea", ssgsea.norm=F)
# identical(ssGSEA[1:100, ], ssGSEA_100)   # T
# 
# ssGSEA[1:3, 1:2]
# BIOCARTA_GRANULOCYTES_PATHWAY BIOCARTA_LYM_PATHWAY
# SIDM00001                      97.38957             116.1215
# SIDM00002                    -348.57552            -206.7634
# SIDM00003                    -377.44541            -222.9439
# 
# ssGSEA_100[1:3, 1:2]
# BIOCARTA_GRANULOCYTES_PATHWAY BIOCARTA_LYM_PATHWAY
# SIDM00001                      97.38957             116.1215
# SIDM00002                    -348.57552            -206.7634
# SIDM00003                    -377.44541            -222.9439


# match(stable_100, stable)   # 56, ..., 937
common_genes = intersect(path_genes, colnames(SANGER_RNA))   # 1506

SANGER_Asinh = SANGER_RNA_Asinh[, common_genes]
SANGER_ZGene = SANGER_RNA[, common_genes] %>% scale %>% as.data.frame
SANGER_ZNorm = SANGER_RNA[, common_genes] %>% t %>% scale %>% t %>% as.data.frame
SANGER_ssGSEA = SANGER_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea")
SANGER_ssGSEA_NN = SANGER_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea", ssgsea.norm=F)
SANGER_GSVA = SANGER_RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva")
SANGER_SING = SANGER_RNA %>% singscore_def(Path_List, filt_genes=T)
SANGER_SING_ST = SANGER_RNA %>% singscore_def(Path_List, filt_genes=T, stableGenes=stable_100)

SANGER_ZNorm_NF = SANGER_RNA %>% t %>% scale %>% t %>% as.data.frame
SANGER_ssGSEA_NF = SANGER_RNA %>% gsva_def(Path_List, filt_genes=F, method="ssgsea")
SANGER_ssGSEA_NN_NF = SANGER_RNA %>% gsva_def(Path_List, filt_genes=F, method="ssgsea", ssgsea.norm=F)
SANGER_GSVA_NF = SANGER_RNA %>% gsva_def(Path_List, filt_genes=F, method="gsva")
SANGER_SING_NF = SANGER_RNA %>% singscore_def(Path_List, filt_genes=F)
SANGER_SING_ST_NF = SANGER_RNA %>% singscore_def(Path_List, filt_genes=F, stableGenes=stable_100)


methods = c("Z_Sample", "ssGSEA", "ssGSEA_Norm_X", 
            "GSVA", "singscore", "stingscore")

# Set the filtering non-pathway genes as the reference data
Corr_Filt_ZNorm = examine_cell_pair(SANGER_ZNorm, SANGER_ZNorm_NF[, common_genes], methods[1])
Corr_Filt_ssGSEA = examine_cell_pair(SANGER_ssGSEA, SANGER_ssGSEA_NF, methods[2])
Corr_Filt_ssGSEA_NN = examine_cell_pair(SANGER_ssGSEA_NN, SANGER_ssGSEA_NN_NF, methods[3])
Corr_Filt_GSVA = examine_cell_pair(SANGER_GSVA, SANGER_GSVA_NF, methods[4])
Corr_Filt_SING = examine_cell_pair(SANGER_SING, SANGER_SING_NF, methods[6])
Corr_Filt_SING_ST = examine_cell_pair(SANGER_SING_ST, SANGER_SING_ST_NF, methods[7])


# Example [SIDM00001]
Ex_Filt_ZNorm = data.frame(Filt_O=as.numeric(SANGER_ZNorm[1, ]), 
                           Filt_X=as.numeric(SANGER_ZNorm_NF[1, ]))
Ex_Filt_ssGSEA = data.frame(Filt_O=as.numeric(SANGER_ssGSEA[1, ]), 
                            Filt_X=as.numeric(SANGER_ssGSEA_NF[1, ]))
Ex_Filt_ssGSEA_NN = data.frame(Filt_O=as.numeric(SANGER_ssGSEA_NN[1, ]), 
                               Filt_X=as.numeric(SANGER_ssGSEA_NN_NF[1, ]))
Ex_Filt_GSVA = data.frame(Filt_O=as.numeric(SANGER_GSVA[1, ]), 
                          Filt_X=as.numeric(SANGER_GSVA_NF[1, ]))
Ex_Filt_SING = data.frame(Filt_O=as.numeric(SANGER_SING[1, ]), 
                          Filt_X=as.numeric(SANGER_SING_NF[1, ]))
Ex_Filt_SING_ST = data.frame(Filt_O=as.numeric(SANGER_SING_ST[1, ]), 
                             Filt_X=as.numeric(SANGER_SING_ST_NF[1, ]))


methods = c("ssGSEA", "GSVA", "singscore")
main = sprintf("%s/Example Scatter [%s, SIDM00001, Genes filtering]", dir, methods)

plot_ex_filt = function(Filt_Ex, main=NULL, save=T, ...) {
  xlab = "SIDM00001 [Pathway Genes]"
  ylab = "SIDM00001 [Total Genes]"
  Filt_Ex %>% plot_def(Filt_O, Filt_X, main=main, xlab=xlab, ylab=ylab, xy_line=T, save=save, ...)
}

Ex_Filt_ZNorm %>% plot_ex_filt(main=main[1])
Ex_Filt_ssGSEA %>% plot_ex_filt(main=main[2])
Ex_Filt_ssGSEA_NN %>% plot_ex_filt(main=main[3])
Ex_Filt_GSVA %>% plot_ex_filt(main=main[4])
Ex_Filt_SING %>% plot_ex_filt(main=main[5])
Ex_Filt_SING_ST %>% plot_ex_filt(main=main[6])



### 2-2. GSVA can calculate the pathway scores differently in different cells? [O]
# Except GSVA, all methods process cells in similar values...

Corr_Diff_ZNorm = SANGER_ZNorm %>% examine_cell(methods[1])
Corr_Diff_ssGSEA = SANGER_ssGSEA %>% examine_cell(methods[2])
Corr_Diff_ssGSEA_NN = SANGER_ssGSEA_NN %>% examine_cell(methods[3])
Corr_Diff_GSVA = SANGER_GSVA %>% examine_cell(methods[4])
Corr_Diff_SING = SANGER_SING %>% examine_cell(methods[5])
Corr_Diff_SING_ST = SANGER_SING_ST %>% examine_cell(methods[6])


Corr_Diff = Reduce(rbind, list(Corr_Diff_ZNorm, Corr_Diff_ssGSEA, Corr_Diff_ssGSEA_NN, 
                               Corr_Diff_GSVA, Corr_Diff_SING, Corr_Diff_SING_ST))

Corr_Diff = Corr_Diff %>% 
  mutate(Method = Method %>% factor(levels=methods))


stat = c("NRMSE", "Corr")
ylab = sprintf("%s [Different cell-pairs]", stat)
dir = mkdir("../../processed_data/cell_data/BIOCARTA/GSVA_Analysis/Cell Distinction")
file = sprintf("%s/Cell-Pair Comparison [%s, SANGER]", dir, stat)

Corr_Diff %>% subset(Cell1!=Cell2) %>% as.data.frame %>% 
  boxplot_def(Method, NRMSE, fill=Method, main=file[1], dpi=1500, 
              ylab=ylab[1], point=F, violin=F, legend=F, save=T, save_svg=T, raster=T)

Corr_Diff %>% subset(Cell1!=Cell2) %>% as.data.frame %>% 
  boxplot_def(Method, Corr, fill=Method, main=file[2], dpi=1500, 
              ylab=ylab[2], point=F, violin=F, legend=F, save=T, save_svg=T, raster=T)


# Example
plot_ex_diff = function(Omics, main=NULL, width=15, height=15, save=T, ...) {
  
  xlab = "SIDM00001 [MEC-1, CLL]"
  ylab = "SIDM00002 [NBsusSR, NB]"
  pcc = Omics %>% t %>% as.data.frame %>% with(cor(SIDM00001, SIDM00002)) %>% round(3)
  sprintf("PCC [SIDM00001-SIDM00002] : %s", pcc) %>% print
  
  Omics %>% t %>% as.data.frame %>% 
    plot_def(SIDM00001, SIDM00002, main=main, alpha=0.5, size=2,
             xlab=xlab, ylab=ylab, width=width, height=height, 
             text_ratio=1.35, xy_line=T, save=save, save_svg=save, raster=T)
}

dir = mkdir("../../processed_data/cell_data/BIOCARTA/GSVA_Analysis/Cell Distinction")
main = sprintf("%s/Example Scatter [%s, SIDM00001 & SIDM00002]", dir, methods)

SANGER_ZNorm %>% plot_ex_diff(main=main[1])
SANGER_ssGSEA %>% plot_ex_diff(main=main[2])
SANGER_ssGSEA_NN %>% plot_ex_diff(main=main[3])
SANGER_GSVA %>% plot_ex_diff(main=main[4])
SANGER_SING %>% plot_ex_diff(main=main[5])
SANGER_SING_ST %>% plot_ex_diff(main=main[6])

file = sprintf("%s/Cell_Distinction.csv", dir)
fwrite(Corr_Diff, file=file)



### 2-3-1. GSVA can alleviate the batch-effects? [O]
# SANGER RNA-Seq vs CCLE RNA-Seq

# all(stable_100 %in% colnames(CCLE_RNA))     # 100 [from 100]
# all(stable_100 %in% colnames(GDSC_RNA))     # 93 [from 100]
all(common_genes %in% colnames(CCLE_RNA))   # 1498 [from 1506]
all(common_genes %in% colnames(GDSC_RNA))   # 1410 [from 1506]
Reduce(intersect, list(stable_100, path_genes, colnames(CCLE_RNA))) %>% length   # 100
Reduce(intersect, list(stable_100, path_genes, colnames(GDSC_RNA))) %>% length   # 93

common_ccle = Reduce(intersect, list(path_genes, colnames(SANGER_RNA), colnames(CCLE_RNA)))
common_gdsc = Reduce(intersect, list(path_genes, colnames(SANGER_RNA), colnames(GDSC_RNA)))

CCLE_Asinh = CCLE_RNA_Asinh[, common_ccle]
CCLE_ZGene = CCLE_RNA[, common_ccle] %>% scale %>% as.data.frame
CCLE_ZNorm = CCLE_RNA[, common_ccle] %>% t %>% scale %>% t %>% as.data.frame
CCLE_ssGSEA = CCLE_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea")
CCLE_ssGSEA_NN = CCLE_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea", ssgsea.norm=F)
CCLE_GSVA = CCLE_RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva")
CCLE_SING = CCLE_RNA %>% singscore_def(Path_List, filt_genes=T)
CCLE_SING_ST = CCLE_RNA %>% singscore_def(Path_List, filt_genes=T, stableGenes=stable_100)

GDSC_ZGene = GDSC_RNA[, common_gdsc] %>% scale %>% as.data.frame
GDSC_ZNorm = GDSC_RNA[, common_gdsc] %>% t %>% scale %>% t %>% as.data.frame
GDSC_ssGSEA = GDSC_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea")
GDSC_ssGSEA_NN = GDSC_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea", ssgsea.norm=F)
GDSC_GSVA = GDSC_RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva")
GDSC_SING = GDSC_RNA %>% singscore_def(Path_List, filt_genes=T)
GDSC_SING_ST = GDSC_RNA %>% singscore_def(Path_List, filt_genes=T, stableGenes=stable_100)


methods = c("No_Norm", "ArcSinh", "Z_Sample", "Z_Gene", "ssGSEA", "ssGSEA_Norm_X", "GSVA", "singscore", "stingscore")

Corr_CCLE_Raw = examine_cell_pair(SANGER_RNA[, common_ccle], CCLE_RNA[, common_ccle], methods[1])
Corr_CCLE_Asinh = examine_cell_pair(SANGER_Asinh[, common_ccle], CCLE_Asinh[, common_ccle], methods[2])
Corr_CCLE_ZNorm = examine_cell_pair(SANGER_ZNorm[, common_ccle], CCLE_ZNorm[, common_ccle], methods[3])
Corr_CCLE_ZGene = examine_cell_pair(SANGER_ZGene[, common_ccle], CCLE_ZGene[, common_ccle], methods[4])
Corr_CCLE_ssGSEA = examine_cell_pair(SANGER_ssGSEA, CCLE_ssGSEA, methods[5])
Corr_CCLE_ssGSEA_NN = examine_cell_pair(SANGER_ssGSEA_NN, CCLE_ssGSEA_NN, methods[6])
Corr_CCLE_GSVA = examine_cell_pair(SANGER_GSVA, CCLE_GSVA, methods[7])
Corr_CCLE_SING = examine_cell_pair(SANGER_SING, CCLE_SING, methods[8])
Corr_CCLE_SING_ST = examine_cell_pair(SANGER_SING_ST, CCLE_SING_ST, methods[9])

Corr_GDSC_Raw = examine_cell_pair(SANGER_RNA[, common_gdsc], GDSC_RNA[, common_gdsc], methods[1])
Corr_GDSC_Asinh = examine_cell_pair(SANGER_Asinh[, common_gdsc], GDSC_RNA[, common_gdsc], methods[2])
Corr_GDSC_ZNorm = examine_cell_pair(SANGER_ZNorm[, common_gdsc], GDSC_ZNorm[, common_gdsc], methods[3])
Corr_GDSC_ZGene = examine_cell_pair(SANGER_ZGene[, common_gdsc], GDSC_ZGene[, common_gdsc], methods[4])
Corr_GDSC_ssGSEA = examine_cell_pair(SANGER_ssGSEA, GDSC_ssGSEA, methods[5])
Corr_GDSC_ssGSEA_NN = examine_cell_pair(SANGER_ssGSEA_NN, GDSC_ssGSEA_NN, methods[6])
Corr_GDSC_GSVA = examine_cell_pair(SANGER_GSVA, GDSC_GSVA, methods[7])
Corr_GDSC_SING = examine_cell_pair(SANGER_SING, GDSC_SING, methods[8])
Corr_GDSC_SING_ST = examine_cell_pair(SANGER_SING_ST, GDSC_SING_ST, methods[9])

col = c("Cell1", "Cell2", "Corr", "Method")
rbind_ = function(df1, df2) rbind(df1[, col], df2[, col])

Corr_SANGER_CCLE = Reduce(rbind_, 
  list(Corr_CCLE_Raw, Corr_CCLE_Asinh, Corr_CCLE_ZNorm, Corr_CCLE_ZGene, 
       Corr_CCLE_ssGSEA, Corr_CCLE_ssGSEA_NN, 
       Corr_CCLE_GSVA, Corr_CCLE_SING, Corr_CCLE_SING_ST)
)   # 14996880

Corr_SANGER_GDSC = Reduce(rbind_, 
  list(Corr_GDSC_Raw, Corr_GDSC_Asinh, Corr_GDSC_ZNorm, Corr_GDSC_ZGene, 
       Corr_GDSC_ssGSEA, Corr_GDSC_ssGSEA_NN, 
       Corr_GDSC_GSVA, Corr_GDSC_SING, Corr_GDSC_SING_ST)
)   # 10077102

colnames(Corr_SANGER_CCLE)[1:2] = c("SANGER", "CCLE")
colnames(Corr_SANGER_GDSC)[1:2] = c("SANGER", "GDSC")

# idx = match(Corr_SANGER_CCLE$CCLE, Anno_Cells$BROAD_ID)
# Anno_Cells$BROAD_ID[idx] %>% is.na %>% sum
# Corr_SANGER_CCLE$CCLE_SANGER_ID = Anno_Cells$BROAD_ID[idx]

Corr_SANGER_CCLE = Corr_SANGER_CCLE %>% 
  mutate(SANGER=as.character(SANGER), CCLE=as.character(CCLE),
         Cell_Pair=ifelse(SANGER==CCLE, "Same", "Different"), 
         Method=factor(Method, level=methods)) %>% 
  mutate(Cell_Pair=Cell_Pair %>% factor(levels=c("Same", "Different")))

Corr_SANGER_GDSC = Corr_SANGER_GDSC %>% 
  mutate(SANGER=as.character(SANGER), GDSC=as.character(GDSC),
         Cell_Pair=ifelse(SANGER==GDSC, "Same", "Different"), 
         Method=factor(Method, level=methods)) %>% 
  mutate(Cell_Pair=Cell_Pair %>% factor(levels=c("Same", "Different")))

Corr_SANGER_CCLE$Cell_Pair %>% table   # Same 9063 [1007*9]
Corr_SANGER_GDSC$Cell_Pair %>% table   # Same 8766 [974*9]



width = 24
height = 16
legend = "Cell-Pair"

dir = mkdir("../../processed_data/cell_data/BIOCARTA/GSVA_Analysis/Batch Correction")
main = c("SANGER TPM & CCLE TPM", "SANGER TPM & GDSC Array")

file1 = sprintf("%s/Cell-Pair Comparison [%s, PCC]", dir, main)
file2 = sprintf("%s/Cell-Pair Comparison [%s, NRMSE]", dir, main)

pos = position_dodge(width=0.85)
ylab = c("Cell-Pair PCC", "Cell-Pair NRMSE [IQR]")
color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 3)]
names(color) = c("Same", "Different")
add = list(scale_fill_manual(values=color))

Corr_SANGER_CCLE %>% as.data.frame %>% 
  boxplot_def(Method, Corr, fill=Cell_Pair, main=file1[1], add=add, pos=pos, 
              ylab=ylab[1], legend=legend, width=width, height=height, 
              alpha=0.9, text_ratio=1.8, hjust=1, vjust=1, dpi=1200, 
              point=F, violin=F, raster=T, save=T, save_svg=T)

Corr_SANGER_GDSC %>% as.data.frame %>% 
  boxplot_def(Method, Corr, fill=Cell_Pair, main=file1[2], add=add, pos=pos, 
              ylab=ylab[1], legend=legend, width=width, height=height, 
              alpha=0.9, text_ratio=1.8, hjust=1, vjust=1, dpi=1200, 
              point=F, violin=F, raster=T, save=T, save_svg=T)

# Corr_SANGER_CCLE %>% as.data.frame %>% 
#   boxplot_def(Method, NRMSE, fill=Cell_Pair, main=file2[2], add=add, pos=pos, 
#               ylab=ylab[2], legend=legend, width=width, height=height, 
#               alpha=0.9, text_ratio=1.65, hjust=1, vjust=1, dpi=1200, 
#               point=F, violin=F, raster=T, save=T, save_svg=T)
# 
# Corr_SANGER_GDSC %>% as.data.frame %>% 
#   boxplot_def(Method, NRMSE, fill=Cell_Pair, main=file2[1], add=add, pos=pos, 
#               ylab=ylab[2], legend=legend, width=width, height=height, 
#               alpha=0.9, text_ratio=1.65, hjust=1, vjust=1, dpi=1200, 
#               point=F, violin=F, raster=T, save=T, save_svg=T)

main = c("SANGER TPM & CCLE TPM", "SANGER TPM & GSDC Array")
file = sprintf("%s/Batch Correction [%s].csv", dir, main)
fwrite(Corr_SANGER_CCLE, file=file[1])
fwrite(Corr_SANGER_GDSC, file=file[2])




### 2-3-2. GSVA can alleviate the batch-effects? [O, Standardization]

# Scaler fitted on SANGER RNA-Seq
Scaler_Raw_S = SANGER_RNA[, common_genes] %>% caret::preProcess()
Scaler_Asinh_S = SANGER_Asinh[, common_genes] %>% caret::preProcess()
Scaler_ZNorm_S = SANGER_ZNorm[, common_genes] %>% caret::preProcess()
Scaler_ZGene_S = SANGER_ZGene[, common_genes] %>% caret::preProcess()

Scaler_Raw_C = SANGER_RNA[, common_ccle] %>% caret::preProcess()
Scaler_Asinh_C = SANGER_Asinh[, common_ccle] %>% caret::preProcess()
Scaler_ZNorm_C = SANGER_ZNorm[, common_ccle] %>% caret::preProcess()
Scaler_ZGene_C = SANGER_ZGene[, common_ccle] %>% caret::preProcess()

Scaler_Raw_G = SANGER_RNA[, common_gdsc] %>% caret::preProcess()
Scaler_Asinh_G = SANGER_Asinh[, common_gdsc] %>% caret::preProcess()
Scaler_ZNorm_G = SANGER_ZNorm[, common_gdsc] %>% caret::preProcess()
Scaler_ZGene_G = SANGER_ZGene[, common_gdsc] %>% caret::preProcess()

Scaler_ssGSEA = SANGER_ssGSEA %>% caret::preProcess()
Scaler_ssGSEA_NN = SANGER_ssGSEA_NN %>% caret::preProcess()
Scaler_GSVA = SANGER_GSVA %>% caret::preProcess()
Scaler_SING = SANGER_SING %>% caret::preProcess()
Scaler_SING_ST = SANGER_SING_ST %>% caret::preProcess()

# Transform SANGER RNA-Seq
scale_list = function(..., Scaler_List=NULL) {
  
  Omics_List = list(...)
  Omics_List_Scaled = list()
  
  if (length(Omics_List)!=length(Scaler_List)) {
    stop("The number of omics and scalers is not identical...")
  }
  
  for (i in 1:length(Omics_List)) {
    Omics_List_Scaled[[i]] = stats::predict(Scaler_List[[i]], Omics_List[[i]])
  }
  return(Omics_List_Scaled)
}

Scaler_List = list(Scaler_Raw_S, Scaler_Asinh_S, Scaler_ZNorm_S, 
                   Scaler_ZGene_S, Scaler_ssGSEA, Scaler_ssGSEA_NN, 
                   Scaler_GSVA, Scaler_SING, Scaler_SING_ST)

Scaler_List_CCLE = list(Scaler_Raw_C, Scaler_Asinh_C, Scaler_ZNorm_C, 
                        Scaler_ZGene_C, Scaler_ssGSEA, Scaler_ssGSEA_NN, 
                        Scaler_GSVA, Scaler_SING, Scaler_SING_ST)

Scaler_List_GDSC = list(Scaler_Raw_G, Scaler_Asinh_G, Scaler_ZNorm_G, 
                        Scaler_ZGene_G, Scaler_ssGSEA, Scaler_ssGSEA_NN, 
                        Scaler_GSVA, Scaler_SING, Scaler_SING_ST)

SANGER_List = scale_list(SANGER_RNA[, common_genes], SANGER_Asinh[, common_genes], 
                         SANGER_ZNorm[, common_genes], SANGER_ZGene[, common_genes], 
                         SANGER_ssGSEA, SANGER_ssGSEA_NN, SANGER_GSVA, 
                         SANGER_SING, SANGER_SING_ST, Scaler_List=Scaler_List)

CCLE_List = scale_list(CCLE_RNA[, common_ccle], CCLE_Asinh[, common_ccle],
                       CCLE_ZNorm[, common_ccle], CCLE_ZGene[, common_ccle], 
                       CCLE_ssGSEA, CCLE_ssGSEA_NN, CCLE_GSVA, 
                       CCLE_SING, CCLE_SING_ST, Scaler_List=Scaler_List_CCLE)

GDSC_List = scale_list(GDSC_RNA[, common_gdsc], GDSC_RNA[, common_gdsc], 
                       GDSC_ZNorm[, common_gdsc], GDSC_ZGene[, common_gdsc], 
                       GDSC_ssGSEA, GDSC_ssGSEA_NN, GDSC_GSVA, 
                       GDSC_SING, GDSC_SING_ST, Scaler_List=Scaler_List_GDSC)

# Recalculate cell-pair correlation
methods_ = c("Raw", "Asinh", "ZNorm", "ZGene", "ssGSEA", "ssGSEA_NN", "GSVA", "SING", "SING_ST")
methods = c("No_Norm", "ArcSinh", "Z_Sample", "Z_Gene", "ssGSEA", "ssGSEA_Norm_X", "GSVA", "singscore", "stingscore")
corr_ccle = sprintf("Corr_CCLE_%s_", methods_)
corr_gdsc = sprintf("Corr_GDSC_%s_", methods_)

for (i in 1:length(methods_)) {
  assign(corr_ccle[i], examine_cell_pair(SANGER_List[[i]], CCLE_List[[i]], methods[i]))
  assign(corr_gdsc[i], examine_cell_pair(SANGER_List[[i]], GDSC_List[[i]], methods[i]))
}

Corr_SANGER_CCLE_Norm = Reduce(rbind, mget(corr_ccle))
Corr_SANGER_GDSC_Norm = Reduce(rbind, mget(corr_gdsc))                                         
colnames(Corr_SANGER_CCLE_Norm)[1:2] = c("SANGER", "CCLE")
colnames(Corr_SANGER_GDSC_Norm)[1:2] = c("SANGER", "GDSC")

# idx = match(Corr_SANGER_CCLE$CCLE, Anno_Cells$BROAD_ID)
# Anno_Cells$BROAD_ID[idx] %>% is.na %>% sum
# Corr_SANGER_CCLE$CCLE_SANGER_ID = Anno_Cells$BROAD_ID[idx]

Corr_SANGER_CCLE_Norm = Corr_SANGER_CCLE_Norm %>% 
  mutate(SANGER=as.character(SANGER), CCLE=as.character(CCLE),
         Cell_Pair=ifelse(SANGER==CCLE, "Same", "Different"), 
         Method=factor(Method, level=methods)) %>% 
  mutate(Cell_Pair=Cell_Pair %>% factor(levels=c("Same", "Different")))

Corr_SANGER_GDSC_Norm = Corr_SANGER_GDSC_Norm %>% 
  mutate(SANGER=as.character(SANGER), GDSC=as.character(GDSC),
         Cell_Pair=ifelse(SANGER==GDSC, "Same", "Different"), 
         Method=factor(Method, level=methods)) %>% 
  mutate(Cell_Pair=Cell_Pair %>% factor(levels=c("Same", "Different")))

Corr_SANGER_CCLE_Norm$Cell_Pair %>% table   # Same 9063 [1007*9]
Corr_SANGER_GDSC_Norm$Cell_Pair %>% table   # Same 8766 [974*9]


dir = mkdir("../../processed_data/cell_data/BIOCARTA/GSVA_Analysis/Batch Correction")
main = c("SANGER TPM & CCLE TPM", "SANGER TPM & GDSC Array")

legend = "Cell-Pair"
file1 = sprintf("%s/Cell-Pair Comparison [%s, PCC (Scaled)]", dir, main)
# file2 = sprintf("%s/Cell-Pair Comparison [%s, NRMSE (Scaled)]", dir, main)

pos = position_dodge(width=0.85)
ylab = c("Cell-Pair PCC", "Cell-Pair NRMSE [IQR]")
color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 3)]
names(color) = c("Same", "Different")
add = list(scale_fill_manual(values=color))

Corr_SANGER_CCLE_Norm %>% as.data.frame %>%
  boxplot_def(Method, Corr, fill=Cell_Pair, main=file1[1], add=add, pos=pos,
              ylab=ylab[1], legend=legend, width=width, height=height,
              alpha=0.9, text_ratio=1.8, hjust=1, vjust=1, dpi=1200,
              point=F, violin=F, raster=T, save=T, save_svg=T)

Corr_SANGER_GDSC_Norm %>% as.data.frame %>%
  boxplot_def(Method, Corr, fill=Cell_Pair, main=file1[2], add=add, pos=pos,
              ylab=ylab[1], legend=legend, width=width, height=height,
              alpha=0.9, text_ratio=1.8, hjust=1, vjust=1, dpi=1200,
              point=F, violin=F, raster=T, save=T, save_svg=T)

# Corr_SANGER_CCLE_Norm %>% as.data.frame %>% 
#   boxplot_def(Method, NRMSE, fill=Cell_Pair, main=file2[2], add=add, pos=pos, 
#               ylab=ylab[2], legend=legend, width=width, height=height, 
#               alpha=0.9, text_ratio=1.65, hjust=1, vjust=1, dpi=1200, 
#               point=F, violin=F, raster=T, save=T, save_svg=T)
# 
# Corr_SANGER_GDSC_Norm %>% as.data.frame %>% 
#   boxplot_def(Method, NRMSE, fill=Cell_Pair, main=file2[1], add=add, pos=pos, 
#               ylab=ylab[2], legend=legend, width=width, height=height, 
#               alpha=0.9, text_ratio=1.65, hjust=1, vjust=1, dpi=1200, 
#               point=F, violin=F, raster=T, save=T, save_svg=T)

main = c("SANGER TPM & CCLE TPM", "SANGER TPM & GSDC Array")
file = sprintf("%s/Batch Correction [%s (Scaled)].csv", dir, main)
fwrite(Corr_SANGER_CCLE_Norm, file=file[1])
fwrite(Corr_SANGER_GDSC_Norm, file=file[2])


pca_chembl = F
if (pca_chembl) {
  # Cell_Drug_Pair
  dir = "../../benchmark_test/_performance_chembl"
  file = sprintf("%s/Cell_Drug_Pair.csv", dir)
  Cell_Drug_Pair = read.csv(file)
  
  # Anno_Cells_ChEMBL
  dir = "../../processed_data/ic50_data/ChEMBL"
  file = sprintf("%s/Anno_Cells_ChEMBL.csv", dir)
  Anno_Cells_ChEMBL = read.csv(file)
  
  cells_chembl = Anno_Cells_ChEMBL %>% 
    subset(Cell_ChEMBL_ID %in% Cell_Drug_Pair$Cell_ChEMBL_ID) %>% 
    pull(SANGER_MODEL_ID) %>% unique   # 404
  
  common_genes_ = Reduce(intersect, list(common_genes, common_ccle, common_gdsc))   # 1408
  RNA_List_Raw = list(SANGER_RNA[, common_genes_], CCLE_RNA[, common_genes_], GDSC_RNA[, common_genes_])
  RNA_List_Asinh = list(SANGER_Asinh[, common_genes_], CCLE_Asinh[, common_genes_], GDSC_RNA[, common_genes_])
  RNA_List_ZNorm = list(SANGER_ZNorm[, common_genes_], CCLE_ZNorm[, common_genes_], GDSC_ZNorm[, common_genes_])
  RNA_List_ZGene = list(SANGER_ZGene[, common_genes_], CCLE_ZGene[, common_genes_], GDSC_ZGene[, common_genes_])
  RNA_List_ssGSEA = list(SANGER_ssGSEA, CCLE_ssGSEA, GDSC_ssGSEA)
  RNA_List_ssGSEA_NN = list(SANGER_ssGSEA_NN, CCLE_ssGSEA_NN, GDSC_ssGSEA_NN)
  RNA_List_GSVA = list(SANGER_GSVA, CCLE_GSVA, GDSC_GSVA)
  RNA_List_SING = list(SANGER_SING, CCLE_SING, GDSC_SING)
  RNA_List_SING_ST = list(SANGER_SING_ST, CCLE_SING_ST, GDSC_SING_ST)
  
  select_chembl = function(df) df[cells_chembl, ]
  RNA_List_Raw = RNA_List_Raw %>% lapply(select_chembl)
  RNA_List_Asinh = RNA_List_Asinh %>% lapply(select_chembl)
  RNA_List_ZNorm = RNA_List_ZNorm %>% lapply(select_chembl)
  RNA_List_ZGene = RNA_List_ZGene %>% lapply(select_chembl)
  RNA_List_ssGSEA = RNA_List_ssGSEA %>% lapply(select_chembl)
  RNA_List_ssGSEA_NN = RNA_List_ssGSEA_NN %>% lapply(select_chembl)
  RNA_List_GSVA = RNA_List_GSVA %>% lapply(select_chembl)
  RNA_List_SING = RNA_List_SING %>% lapply(select_chembl)
  RNA_List_SING_ST = RNA_List_SING_ST %>% lapply(select_chembl)
  
  width = 18.6
  height = 13.5
  db_list = c("SANGER", "CCLE", "GDSC")
  
  methods = c("No_Norm", "ArcSinh", "Z_Sample", "Z_Gene", "ssGSEA", "ssGSEA_Norm_X", "GSVA", "singscore", "stingscore")
  dir = mkdir("../../benchmark_test/_performance_chembl/PC2 Plot [SANGER & CCLE & GDSC]")
  main = sprintf("%s/PC2 Plot of SANGER & CCLE & GDSC [%s]", dir, methods)
  
  shape = c(21, 22, 24)
  color = c("brown1", "seagreen3", "royalblue1")
  
  add = list(scale_color_manual(values=color), 
             scale_shape_manual(values=shape), 
             theme(legend.key.size=unit(1, "cm")), 
             guides(shape=guide_legend(override.aes=list(size=3, alpha=1))))
  
  PCA_Raw = RNA_List_Raw %>% plot_pca_batch(
    db_list, main=main[1], size=1.5, alpha=0.5, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, height=height, add=add, save=T)
  PCA_Asinh = RNA_List_Asinh %>% plot_pca_batch(
    db_list, main=main[2], size=1.5, alpha=0.5, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, height=height, add=add, save=T)
  PCA_ZNorm = RNA_List_ZNorm %>% plot_pca_batch(
    db_list, main=main[3], size=1.5, alpha=0.5, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, height=height, add=add, save=T)
  PCA_ZGene = RNA_List_ZGene %>% plot_pca_batch(
    db_list, main=main[4], size=1.5, alpha=0.5, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, height=height, add=add, save=T)
  PCA_ssGSEA = RNA_List_ssGSEA %>% plot_pca_batch(
    db_list, main=main[5], size=1.5, alpha=0.5, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, height=height, add=add, save=T)
  PCA_ssGSEA_NN = RNA_List_ssGSEA_NN %>% plot_pca_batch(
    db_list, main=main[6], size=1.5, alpha=0.5, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, height=height, add=add, save=T)
  PCA_GSVA = RNA_List_GSVA %>% plot_pca_batch(
    db_list, main=main[7], size=1.5, alpha=0.5, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, height=height, add=add, save=T)
  PCA_SING = RNA_List_SING %>% plot_pca_batch(
    db_list, main=main[8], size=1.5, alpha=0.5, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, height=height, add=add, save=T)
  PCA_SING_ST = RNA_List_SING_ST %>% plot_pca_batch(
    db_list, main=main[9], size=1.5, alpha=0.5, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=20, width=width, height=height, add=add, save=T)
  
  
  width = 24
  height = 16
  legend = "Cell-Pair"
  
  dir = "../../processed_data/cell_data/BIOCARTA/GSVA_Analysis/Batch Correction"
  main = c("SANGER TPM & CCLE TPM", "SANGER TPM & GDSC Array")
  file1 = sprintf("%s/Cell-Pair Comparison [%s for ChEMBL, PCC]", dir, main)
  file2 = sprintf("%s/Cell-Pair Comparison [%s for ChEMBL, PCC (Scaled)]", dir, main)
  
  pos = position_dodge(width=0.85)
  ylab = c("Cell-Pair PCC", "Cell-Pair NRMSE [IQR]")
  color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 3)]
  names(color) = c("Same", "Different")
  add = list(scale_fill_manual(values=color))
  
  Corr_SANGER_CCLE_Norm = Corr_SANGER_CCLE_Norm %>% as.data.table
  Corr_SANGER_GDSC_Norm = Corr_SANGER_GDSC_Norm %>% as.data.table
  
  Corr_List = list()
  Corr_List[[1]] = Corr_SANGER_CCLE[SANGER %in% cells_chembl & CCLE %in% cells_chembl]
  Corr_List[[2]] = Corr_SANGER_GDSC[SANGER %in% cells_chembl & GDSC %in% cells_chembl]
  Corr_List[[3]] = Corr_SANGER_CCLE_Norm[SANGER %in% cells_chembl & CCLE %in% cells_chembl]
  Corr_List[[4]] = Corr_SANGER_GDSC_Norm[SANGER %in% cells_chembl & GDSC %in% cells_chembl]
  
  Corr_List[[1]] %>% as.data.frame %>% 
    boxplot_def(Method, Corr, fill=Cell_Pair, main=file1[1], add=add, pos=pos, 
                ylab=ylab[1], legend=legend, width=width, height=height, 
                alpha=0.9, text_ratio=1.8, hjust=1, vjust=1, dpi=1200, 
                point=T, violin=F, raster=T, save=T, save_svg=T)
  
  Corr_List[[2]] %>% as.data.frame %>% 
    boxplot_def(Method, Corr, fill=Cell_Pair, main=file1[2], add=add, pos=pos, 
                ylab=ylab[1], legend=legend, width=width, height=height, 
                alpha=0.9, text_ratio=1.8, hjust=1, vjust=1, dpi=1200, 
                point=T, violin=F, raster=T, save=T, save_svg=T)
  
  Corr_List[[3]] %>% as.data.frame %>%
    boxplot_def(Method, Corr, fill=Cell_Pair, main=file2[1], add=add, pos=pos,
                ylab=ylab[1], legend=legend, width=width, height=height,
                alpha=0.9, text_ratio=1.8, hjust=1, vjust=1, dpi=1200,
                point=T, violin=F, raster=T, save=T, save_svg=T)
  
  Corr_List[[4]] %>% as.data.frame %>%
    boxplot_def(Method, Corr, fill=Cell_Pair, main=file2[2], add=add, pos=pos,
                ylab=ylab[1], legend=legend, width=width, height=height,
                alpha=0.9, text_ratio=1.8, hjust=1, vjust=1, dpi=1200,
                point=T, violin=F, raster=T, save=T, save_svg=T)
  
  col_old1 = c("SANGER", "CCLE", "Corr", "Method", "Cell_Pair")
  col_old2 = c("SANGER", "GDSC", "Corr", "Method", "Cell_Pair")
  col_new = c("Cell_SANGER", "Cell_Other_DB", "Corr", "Method", "Cell_Pair")
  
  setnames(Corr_List[[1]], col_old1, col_new)
  setnames(Corr_List[[2]], col_old2, col_new)
  setnames(Corr_List[[3]], col_old1, col_new)
  setnames(Corr_List[[4]], col_old2, col_new)
  
  Corr_List[[1]]$Other_DB = "CCLE"
  Corr_List[[2]]$Other_DB = "GDSC"
  Corr_List[[3]]$Other_DB = "CCLE_After_Norm"
  Corr_List[[4]]$Other_DB = "GDSC_After_Norm"
  
  Corr_List = rbindlist(Corr_List)
  setkey(Corr_List, Cell_SANGER, Cell_Other_DB, Method, Cell_Pair, Other_DB)
  
  Corr_List = dcast(
    Corr_List,
    Cell_SANGER+Cell_Other_DB+Method+Cell_Pair~Other_DB,
    value.var = "Corr"
  )
}




### 2-4. Complete Pathway-Pathway correlation network [RNA, k=5]

RNA_Corr = SANGER_GSVA %>% cor %>% reshape2::melt() %>% as.data.frame
colnames(RNA_Corr) = c("Pathway1", "Pathway2", "Corr")
RNA_Corr = RNA_Corr %>% subset(Pathway1!=Pathway2)   # 36672 [192 x 191]

RNA_Corr_KNN3 = RNA_Corr %>% knn_graph(col_stat="Corr", mode="max", k=3)   # 576
RNA_Corr_KNN5 = RNA_Corr %>% knn_graph(col_stat="Corr", mode="max", k=5)   # 960
RNA_Corr_KNN7 = RNA_Corr %>% knn_graph(col_stat="Corr", mode="max", k=7)   # 1344

# RNA_Corr$Corr %>% hist
# RNA_Corr_KNN3$Corr %>% hist
# RNA_Corr_KNN5$Corr %>% hist


# Save files
main = c("SANGER", "GDSC", "CCLE")
dir = mkdir("../../processed_data/cell_data/BIOCARTA")

file = sprintf("%s/%s_RNA_GSVA.csv", dir, main)
write.csv(SANGER_GSVA, file=file[1], row.names=T)
write.csv(GDSC_RNA_GSVA, file=file[2], row.names=T)
write.csv(CCLE_RNA_GSVA, file=file[3], row.names=T)

file = sprintf("%s/%s_RNA_ssGSEA.csv", dir, main)
write.csv(SANGER_ssGSEA, file=file[1], row.names=T)
write.csv(GDSC_RNA_ssGSEA, file=file[2], row.names=T)
write.csv(CCLE_RNA_ssGSEA, file=file[3], row.names=T)

file = sprintf("%s/%s_RNA_SING.csv", dir, main)
write.csv(SANGER_SING, file=file[1], row.names=T)
write.csv(GDSC_RNA_SING, file=file[2], row.names=T)
write.csv(CCLE_RNA_SING, file=file[3], row.names=T)

dir = mkdir("../../processed_data/net_data/BIOCARTA")
file = sprintf("%s/RNA_Corr_KNN%s.csv", dir, c(3, 5, 7))
write.csv(RNA_Corr_KNN3, file=file[1], row.names=F)
write.csv(RNA_Corr_KNN5, file=file[2], row.names=F)
write.csv(RNA_Corr_KNN7, file=file[3], row.names=F)

file = sprintf("%s/RNA_Corr.csv", dir)
write.csv(RNA_Corr, file=file, row.names=F)



### 2-5. Compress STRING & RegNetwork into Pathway-Pathway Association network [k=5]
# [2015] Uncovering disease-disease relationships through the incomplete interactome

dir = "../../processed_data/net_data/STRING"
file= sprintf("%s/STRING_Filt_Ent.csv", dir)
STRING_700 = read.csv(file)

dir = "../../processed_data/net_data/RegNetwork"
file= sprintf("%s/RegNetwork_Ent.csv", dir)
RegNetwork = read.csv(file)

STRING_700 = STRING_700 %>% mutate(Node1=as.character(Node1), 
                                   Node2=as.character(Node2))
RegNetwork = RegNetwork %>% mutate(Node1=as.character(Node1), 
                                   Node2=as.character(Node2))
STRING_900 = STRING_700 %>% subset(Combined_Score>=900)   # 479908 > 236392

col = c("Node1", "Node2")
row_overlap(STRING_700, RegNetwork, col=col)   # 10036 / 371910 [2.70%]

STRING_700_SS = STRING_700 %>% separation_score(Path_List)
# Pathway size range : [10, 351]
# Total genes utilized : 92.06% [4939/5365]

STRING_900_SS = STRING_900 %>% separation_score(Path_List)
# Pathway size range : [9, 345]
# Total Genes utilized : 85.83% [4605/5365]

RegNetwork_SS = RegNetwork %>% separation_score(Path_List)
# Pathway size range : [10, 349]
# Total Genes utilized : 94.28% [5058/5365]


# [Confirmation] Do separation scores truly represent the distance of pathway-pathway?
net_type = c("STRING_700", "STRING_900", "RegNetwork")
dir = mkdir("../../processed_data/net_data/BIOCARTA/Network_Analysis")
main = sprintf("%s/Separation Score & Overlap Ratio [%s]", dir, net_type)

xlab = "Separation Score"
ylab = "Overlap Ratio"

# STRING_700_SS %>% 
#   plot_def(Dist, Ratio_Overlap, main=main[1], xlab=xlab, ylab=ylab, 
#            alpha=0.5, axis_tl=20, save=T)   # [△]
# STRING_900_SS %>% 
#   plot_def(Dist, Ratio_Overlap, main=main[2], xlab=xlab, ylab=ylab, 
#            alpha=0.5, axis_tl=20, save=T)   # [△]
# RegNetwork_SS %>% 
#   plot_def(Dist, Ratio_Overlap, main=main[3], xlab=xlab, ylab=ylab, 
#            alpha=0.5,axis_tl=20, save=T)   # [O]

col = list(fill="lightgray")
margin = margin(10, 10, 10, 10, unit="pt")
font_label = font("xylab", size=30, margin=margin)
font_text = font("xy.text", size=22.5, color="grey30", margin=margin)

font_ = font_label + font_text
main = sprintf("%s/Separation Score & Overlap Ratio [%s, with histogram]", dir, net_type)

pl1 = STRING_700_SS %>% 
  ggscatterhist("Dist", "Ratio_Overlap", xlab=xlab, ylab=ylab, 
                alpha=0.25, size=1.5, margin.params=col)
pl2 = STRING_900_SS %>% 
  ggscatterhist("Dist", "Ratio_Overlap", xlab=xlab, ylab=ylab, 
                alpha=0.25, size=1.5, margin.params=col)
pl3 = RegNetwork_SS %>% 
  ggscatterhist("Dist", "Ratio_Overlap", xlab=xlab, ylab=ylab, 
                alpha=0.25, size=1.5, margin.params=col)

pl1$sp = pl1$sp + font_
pl2$sp = pl2$sp + font_
pl3$sp = pl3$sp + font_

suppressMessages(library(ggrastr))
pl1$sp$layers[[1]] = pl1$sp$layers[[1]] %>% rasterise(dpi=1200, dev="ragg_png")
pl2$sp$layers[[1]] = pl2$sp$layers[[1]] %>% rasterise(dpi=1200, dev="ragg_png")
pl3$sp$layers[[1]] = pl3$sp$layers[[1]] %>% rasterise(dpi=1200, dev="ragg_png")

pl1 %>% save_fig_ggpubr(main=main[1], width=15, height=15, dpi=1200, svg=T)
pl2 %>% save_fig_ggpubr(main=main[2], width=15, height=15, dpi=1200, svg=T)
pl3 %>% save_fig_ggpubr(main=main[3], width=15, height=15, dpi=1200, svg=T)

# Some pathways in RegNetwork have negative separative scores with no gene-overlap ratios...
# They 
STRING_700_SS$Dist %>% na.omit %>% range %>% round(3)   # 0.016   2.711
STRING_900_SS$Dist %>% na.omit %>% range %>% round(3)   # 0.020   4.337
RegNetwork_SS$Dist %>% na.omit %>% range %>% round(3)   # -0.827  1.972

RegNetwork_SS %>% subset(Dist<0) %>% 
  plot_def(Dist, Dist_P1_P2, color=Ratio_Overlap, size=Ratio_Overlap, alpha=0.5, save=F)

RegNetwork_SS %>% subset(Dist<0) %>% 
  plot_def(Dist, Ratio_Overlap, color=Dist_P1_P2, size=Dist_P1_P2, alpha=0.25, save=F)

RegNetwork_SS %>% subset(Dist<0) %>% 
  plot_def(Dist_P1_P2, Ratio_Overlap, color=Dist, size=Dist, alpha=0.25, save=F)

RegNetwork_SS %>% subset(Dist<0) %>% 
  plot_def(Dist_P1_P2, Dist_P1+Dist_P2, color=Dist, size=Ratio_Overlap, alpha=0.25, save=F)


# KNN Graphs
STRING_700_KNN3 = STRING_700_SS %>% knn_graph_ovl(k=3)   # 576 [192*3]
STRING_900_KNN3 = STRING_900_SS %>% knn_graph_ovl(k=3)   # 576 [192*3]
RegNetwork_KNN3 = RegNetwork_SS %>% knn_graph_ovl(k=3)   # 568 [192*3-]

STRING_700_KNN5 = STRING_700_SS %>% knn_graph_ovl(k=5)   # 960 [192*5]
STRING_900_KNN5 = STRING_900_SS %>% knn_graph_ovl(k=5)   # 960 [192*5]
RegNetwork_KNN5 = RegNetwork_SS %>% knn_graph_ovl(k=5)   # 933 [192*5-]

STRING_700_KNN7 = STRING_700_SS %>% knn_graph_ovl(k=7)   # 1344 [192*7]
STRING_900_KNN7 = STRING_900_SS %>% knn_graph_ovl(k=7)   # 1344 [192*7]
RegNetwork_KNN7 = RegNetwork_SS %>% knn_graph_ovl(k=7)   # 1289 [192*7-]

# Cf. KNN Graphs with Maximum Separation Scores [Ablation Test]
STRING_700_KNN5_Inv = STRING_700_SS %>% knn_graph_ovl(k=5, inverse=T)
STRING_900_KNN5_Inv = STRING_900_SS %>% knn_graph_ovl(k=5, inverse=T)
RegNetwork_KNN5_Inv = RegNetwork_SS %>% knn_graph_ovl(k=5, inverse=T)


# KNN Graphs [Random]
seed = 2021
Random_KNN3 = knn_graph_rand(Path_List, k=3, seed=seed)
Random_KNN5 = knn_graph_rand(Path_List, k=5, seed=seed)
Random_KNN7 = knn_graph_rand(Path_List, k=7, seed=seed)

# Random_KNN3_ = knn_graph_rand(Path_List, k=3, seed=seed)
# Random_KNN5_ = knn_graph_rand(Path_List, k=5, seed=seed)
# Random_KNN7_ = knn_graph_rand(Path_List, k=7, seed=seed)
# 
# identical(Random_KNN3, Random_KNN3_)   # T
# identical(Random_KNN5, Random_KNN5_)   # T
# identical(Random_KNN7, Random_KNN7_)   # T

net_names = c("STRING_700", "STRING_900", "RegNetwork", "RNA_Corr", "Random")
PPA_Info = examine_graph_list(STRING_700_KNN5, STRING_900_KNN5, RegNetwork_KNN5, 
                              RNA_Corr_KNN5, Random_KNN5, net_names=net_names)


# All graphs have no multiple subgraphs
suppressMessages(library(ggnetwork))
main = c("GSVA_Corr", "RegNetwork", "STRING_700", "STRING_900")
dir = mkdir("../../processed_data/net_data/BIOCARTA/Network_Analysis")
file = sprintf("%s/KNN5 Graph [%s]", dir, main)

label_corr = RNA_Corr_KNN5 %>% choose_freq_path(top=20)
label_regnet = RegNetwork_KNN5 %>% choose_freq_path(top=20)
label_str_700 = STRING_700_KNN5 %>% choose_freq_path(top=20)
label_str_900 = STRING_900_KNN5 %>% choose_freq_path(top=20)

color = "royalblue1"
add = list(geom_nodetext_repel(aes(label=name), max.overlaps=8, force=2, size=3.2))
RNA_Corr_KNN5 %>% ggnet2_def(main=file[1], color=color, label=label_corr, add=add, dpi=1500, save=T)
RegNetwork_KNN5 %>% ggnet2_def(main=file[2], color=color, label=label_regnet, add=add, dpi=1500, save=T)
STRING_700_KNN5 %>% ggnet2_def(main=file[3], color=color, label=label_str_700, add=add, dpi=1500, save=T)
STRING_900_KNN5 %>% ggnet2_def(main=file[4], color=color, label=label_str_900, add=add, dpi=1500, save=T)

# Similarity between KNN Networks [a little]
row_overlap(STRING_700_KNN5, STRING_900_KNN5)   # 694 / 960 [72.29%]
row_overlap(STRING_700_KNN5, RegNetwork_KNN5)   # 226 / 933 [24.22%]
row_overlap(STRING_700_KNN5, RNA_Corr_KNN5)     # 400 / 960 [41.67%]

net_names = c("Random", "GSVA_Corr", "RegNetwork", "STRING_700", "STRING_900")
KNN3_List = list(Random_KNN3, RNA_Corr_KNN3, RegNetwork_KNN3, STRING_700_KNN3, STRING_900_KNN3)
KNN5_List = list(Random_KNN5, RNA_Corr_KNN5, RegNetwork_KNN5, STRING_700_KNN5, STRING_900_KNN5)
KNN7_List = list(Random_KNN7, RNA_Corr_KNN7, RegNetwork_KNN7, STRING_700_KNN7, STRING_900_KNN7)

knn = c("KNN3", "KNN5", "KNN7")
dir = "../../processed_data/net_data/BIOCARTA/Network_Analysis"
main = sprintf("%s/KNN Graph Similarity [%s]", dir, knn)

KNN3_Sim = row_overlap_heatmap(KNN3_List, net_names=net_names, main=main[1], save=T)
KNN5_Sim = row_overlap_heatmap(KNN5_List, net_names=net_names, main=main[2], save=T)
KNN7_Sim = row_overlap_heatmap(KNN7_List, net_names=net_names, main=main[3], save=T)

file = sprintf("%s/PPA_Information.csv", dir)
write.csv(PPA_Info, file=file, row.names=F)




##### 4. Save Files

dir = "../../processed_data/net_data/BIOCARTA"

main = c("STRING_700", "STRING_900", "RegNetwork")
file = sprintf("%s/SS_%s.csv", dir, main)
write.csv(STRING_700_SS, file=file[1], row.names=F)
write.csv(STRING_900_SS, file=file[2], row.names=F)
write.csv(RegNetwork_SS, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_STRING_700.csv", dir, c(3, 5, 7))
write.csv(STRING_700_KNN3, file=file[1], row.names=F)
write.csv(STRING_700_KNN5, file=file[2], row.names=F)
write.csv(STRING_700_KNN7, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_STRING_900.csv", dir, c(3, 5, 7))
write.csv(STRING_900_KNN3, file=file[1], row.names=F)
write.csv(STRING_900_KNN5, file=file[2], row.names=F)
write.csv(STRING_900_KNN7, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_RegNetwork.csv", dir, c(3, 5, 7))
write.csv(RegNetwork_KNN3, file=file[1], row.names=F)
write.csv(RegNetwork_KNN5, file=file[2], row.names=F)
write.csv(RegNetwork_KNN7, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_Random.csv", dir, c(3, 5, 7))
write.csv(Random_KNN3, file=file[1], row.names=F)
write.csv(Random_KNN5, file=file[2], row.names=F)
write.csv(Random_KNN7, file=file[3], row.names=F)

main = c("STRING_700", "STRING_900", "RegNetwork")
file = sprintf("%s/KNN5_%s_Inv.csv", dir, main)
write.csv(STRING_700_KNN5_Inv, file=file[1], row.names=F)
write.csv(STRING_900_KNN5_Inv, file=file[2], row.names=F)
write.csv(RegNetwork_KNN5_Inv, file=file[3], row.names=F)

col = c("Pathway1", "Pathway2")
RNA_Corr_KNN5_Inv = RNA_Corr %>% knn_graph(col_stat="Corr", mode="min", k=5)   # 960
Edge_Inv = Reduce(rbind, list(STRING_900_KNN5_Inv[, col] %>% mutate(Edge_Type="STRING_900"), 
                              RegNetwork_KNN5_Inv[, col] %>% mutate(Edge_Type="RegNetwork"), 
                              RNA_Corr_KNN5_Inv[, col] %>% mutate(Edge_Type="RNA_Corr")))

file = sprintf("%s/KNN5_RNA_Corr_Inv.csv", dir)
write.csv(RNA_Corr_KNN5_Inv, file=file, row.names=F)

file = sprintf("%s/KNN5_STR9_Reg_Corr_Inv.csv", dir)
write.csv(Edge_Inv, file=file, row.names=F)


file = sprintf("%s/Pathway_Analysis.RData", dir)
save(RegNetwork_SS, STRING_700_SS, STRING_900_SS, RNA_Corr,
     Random_KNN3, Random_KNN5, Random_KNN7,
     STRING_700_KNN3, STRING_900_KNN3, RegNetwork_KNN3, RNA_Corr_KNN3, 
     STRING_700_KNN5, STRING_900_KNN5, RegNetwork_KNN5, RNA_Corr_KNN5, 
     STRING_700_KNN7, STRING_900_KNN7, RegNetwork_KNN7, RNA_Corr_KNN7, file=file)


# cf. Subcategories in STRING PPI network [Physical vs Non-Physical]
dir = "../../raw_data/STRING"
file = sprintf("%s/9606.protein.physical.links.v11.5.txt", dir)
STRING_Phy = fread(file)   # 1991832

dir = "../../processed_data/net_data/STRING"
file = sprintf("%s/Anno_Genes.csv", dir)
STRING_Anno = read.csv(file)

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)

colnames(STRING_Phy) = c("Node1", "Node2", "Combined_Score")
STRING_Phy = STRING_Phy %>% subset(Combined_Score>=900)   # 1991832 > 83896

# Ensembl ID > Symbol
idx1 = match(STRING_Phy$Node1, STRING_Anno$ENSEMBL_ID)
idx2 = match(STRING_Phy$Node2, STRING_Anno$ENSEMBL_ID)

STRING_Phy_Sym = STRING_Phy %>% 
  mutate(Node1=STRING_Anno$SYMBOL[idx1],
         Node2=STRING_Anno$SYMBOL[idx2])

STRING_Phy_Sym = STRING_Phy_Sym[complete.cases(STRING_Phy_Sym), ]   # 83896 > 83896

# Symbol > Entrez ID
idx1 = match(STRING_Phy_Sym$Node1, Anno_Genes$HGNC_SYMBOL)
idx2 = match(STRING_Phy_Sym$Node2, Anno_Genes$HGNC_SYMBOL)

STRING_Phy_Ent = STRING_Phy_Sym %>% 
  mutate(Node1=Anno_Genes$ENTREZ_ID[idx1] %>% as.character,
         Node2=Anno_Genes$ENTREZ_ID[idx2] %>% as.character)

STRING_Phy_Ent = STRING_Phy_Ent[complete.cases(STRING_Phy_Ent), ]   # 83896 > 80092

STRING_Phy_Ent$Node1 = STRING_Phy_Ent$Node1 %>% as.numeric
STRING_Phy_Ent$Node2 = STRING_Phy_Ent$Node2 %>% as.numeric
STRING_Rest_Ent = setdiff(STRING_900[, 1:2], STRING_Phy_Ent[, 1:2])   # 236392 > 156678 [-79714]

# STRING with physical interaction
net_genes = unique(unlist(STRING_Phy_Ent$Node1, STRING_Phy_Ent$Node2))
sum(unique(unlist(Path_List)) %in% net_genes)   # 1255 / 1509
sapply(Path_List, function(x) length(intersect(x, net_genes))) %>% range   # [0, 64]

STRING_Phy_SS = STRING_Phy_Ent %>% separation_score(Path_List)
STRING_Phy_KNN5 = STRING_Phy_SS %>% knn_graph_ovl(k=5)

# STRING without physical interaction
net_genes = unique(unlist(STRING_Rest_Ent$Node1, STRING_Rest_Ent$Node2))
sum(unique(unlist(Path_List)) %in% net_genes)   # 1449 / 1509
sapply(Path_List, function(x) length(intersect(x, net_genes))) %>% range   # [4, 81]

STRING_Rest_SS = STRING_Rest_Ent %>% separation_score(Path_List)
STRING_Rest_KNN5 = STRING_Rest_SS %>% knn_graph_ovl(k=5)

dir = "../../processed_data/net_data/BIOCARTA"
file = sprintf("%s/KNN5_STRING_900_%s.csv", dir, c("Physical", "Rest"))
write.csv(STRING_Phy_KNN5, file=file[1], row.names=F)
write.csv(STRING_Rest_KNN5, file=file[2], row.names=F)

col = c("Pathway1", "Pathway2")
Edge_KNN5_Phy = Reduce(rbind, list(STRING_Phy_KNN5[, col] %>% mutate(Edge_Type="STRING_900"), 
                                   RegNetwork_KNN5[, col] %>% mutate(Edge_Type="RegNetwork"), 
                                   RNA_Corr_KNN5[, col] %>% mutate(Edge_Type="RNA_Corr")))

Edge_KNN5_Rest = Reduce(rbind, list(STRING_Rest_KNN5[, col] %>% mutate(Edge_Type="STRING_900"), 
                                    RegNetwork_KNN5[, col] %>% mutate(Edge_Type="RegNetwork"), 
                                    RNA_Corr_KNN5[, col] %>% mutate(Edge_Type="RNA_Corr")))

dir = "../../processed_data/net_data/BIOCARTA"
file1 = sprintf("%s/KNN5_STR9P_Reg_Corr.csv", dir)
file2 = sprintf("%s/KNN5_STR9R_Reg_Corr.csv", dir)
write.csv(Edge_KNN5_Phy, file=file1, row.names=F)
write.csv(Edge_KNN5_Rest, file=file2, row.names=F)




# cf. Add PCN Weights
softmax = function(x, temp=1, na_value=0) {
  if (all(is.na(x))) {
    x = rep(1/length(x), length(x))
  } else {
    x = exp(x / temp)
    x = x / sum(x, na.rm=T)
    x[is.na(x)] = na_value
  }
  return(x)
}

# RNA_Corr_KNN5 = RNA_Corr_KNN5 %>% group_by(Pathway1) %>%
#   mutate(Weight=softmax(Corr)) %>% as.data.frame
# RegNetwork_KNN5 = RegNetwork_KNN5 %>% group_by(Pathway1) %>%
#   mutate(Weight=softmax(-Dist)) %>% as.data.frame
# STRING_700_KNN5 = STRING_700_KNN5 %>% group_by(Pathway1) %>%
#   mutate(Weight=softmax(-Dist)) %>% as.data.frame
# STRING_900_KNN5 = STRING_900_KNN5 %>% group_by(Pathway1) %>%
#   mutate(Weight=softmax(-Dist)) %>% as.data.frame
# 
# STRING_900_SS = STRING_900_SS %>% group_by(Pathway1) %>%
#   mutate(Weight=softmax(-Dist)) %>% as.data.frame
# 
# net_names = c("GSVA_Corr", "RegNetwork", "STRING_700", "STRING_900")
# main = sprintf("Histogram of Edge Weight [%s]", net_names)
# 
# RNA_Corr_KNN5$Weight %>% hist_def(main=main[1], save=T)
# RegNetwork_KNN5$Weight %>% hist_def(main=main[2], save=T)
# STRING_700_KNN5$Weight %>% hist_def(main=main[3], save=T)
# STRING_900_KNN5$Weight %>% hist_def(main=main[4], save=T)
# 
# col = c("Pathway1", "Pathway2", "Weight")
# Edge_KNN5W = Reduce(rbind, list(STRING_900_KNN5[, col] %>% mutate(Edge_Type="STRING_900"),
#                                 RegNetwork_KNN5[, col] %>% mutate(Edge_Type="RegNetwork"),
#                                 RNA_Corr_KNN5[, col] %>% mutate(Edge_Type="RNA_Corr")))


col1 = c("Pathway1", "Pathway2", "Dist")
col2 = c("Pathway1", "Pathway2", "Dist", "Ratio_Overlap")
col3 = c("Pathway1", "Pathway2", "Corr")

by = c("Pathway1"="Pathway1", "Pathway2"="Pathway2")
merge_ = function(df1, df2) full_join(df1, df2, by=by)
Edge_Feat = Reduce(merge_, list(STRING_900_SS[, col1], RegNetwork_SS[, col2], RNA_Corr[, col3]))

Edge_Feat = Edge_Feat %>% 
  mutate(Weight_STR9=-Dist.x, Weight_Reg=-Dist.y, 
         Weight_Corr=Corr, Weight_OVL=Ratio_Overlap) %>% 
  subset(select=-c(Dist.x, Dist.y, Corr, Ratio_Overlap)) %>% as.data.frame

Edge_Feat$Weight_STR9 %>% is.na %>% sum   # 0
Edge_Feat$Weight_Reg %>% is.na %>% sum    # 30450
Edge_Feat$Weight_Corr %>% is.na %>% sum   # 0
Edge_Feat$Weight_OVL %>% is.na %>% sum    # 0

Edge_Feat$Weight_STR9 %>% range(na.rm=T) %>% round(3)   # -3.164 -0.006
Edge_Feat$Weight_Reg %>% range(na.rm=T) %>% round(3)    # -3.017  0.342
Edge_Feat$Weight_Corr %>% range(na.rm=T) %>% round(3)   # -0.672  0.969
Edge_Feat$Weight_OVL %>% range(na.rm=T) %>% round(3)    # 0.000 0.846

Edge_Feat$Weight_Reg[is.na(Edge_Feat$Weight_Reg)] = min(Edge_Feat$Weight_Reg, na.rm=T)

col = c("Pathway1", "Pathway2")
Edge_KNN5 = Reduce(rbind, list(STRING_900_KNN5[, col] %>% mutate(Edge_Type="STRING_900"),
                               RegNetwork_KNN5[, col] %>% mutate(Edge_Type="RegNetwork"),
                               RNA_Corr_KNN5[, col] %>% mutate(Edge_Type="RNA_Corr")))

by = c("Pathway1"="Pathway1", "Pathway2"="Pathway2")
Edge_KNN5W_Plus = left_join(Edge_KNN5, Edge_Feat, by=by)

Edge_KNN5W_Plus$Weight_STR9 %>% range(na.rm=T) %>% round(3)   # -2.715 -0.006
Edge_KNN5W_Plus$Weight_Reg %>% range(na.rm=T) %>% round(3)    # -3.017  0.342
Edge_KNN5W_Plus$Weight_Corr %>% range(na.rm=T) %>% round(3)   # -0.634  0.969
Edge_KNN5W_Plus$Weight_OVL %>% range(na.rm=T) %>% round(3)    # 0.000 0.846

Edge_KNN5W_Plus$Weight_STR9 %>% hist_def(main="W_STR9", save=T)
Edge_KNN5W_Plus$Weight_Reg %>% hist_def(main="W_Reg", save=T)
Edge_KNN5W_Plus$Weight_Corr %>% hist_def(main="W_Corr", save=T)
Edge_KNN5W_Plus$Weight_OVL %>% hist_def(main="W_OVL", save=T)


robust_scaler = function(x) (x-median(x)) / IQR(x)

Edge_Feat_ = Edge_Feat
Edge_Feat_$Weight_STR9 = Edge_Feat_$Weight_STR9 %>% robust_scaler
Edge_Feat_$Weight_Reg = Edge_Feat_$Weight_Reg %>% robust_scaler
Edge_KNN5W_Plus_ = left_join(Edge_KNN5, Edge_Feat_, by=by)

Edge_KNN5W_Plus_$Weight_STR9 %>% range(na.rm=T) %>% round(3)   # -2.933  1.256
Edge_KNN5W_Plus_$Weight_Reg %>% range(na.rm=T) %>% round(3)    # -0.783  0.533
Edge_KNN5W_Plus_$Weight_Corr %>% range(na.rm=T) %>% round(3)   # -0.634  0.969
Edge_KNN5W_Plus_$Weight_OVL %>% range(na.rm=T) %>% round(3)    # 0.000 0.846

Edge_KNN5W_Plus$Weight_STR9 %>% hist_def(main="W_STR9", save=T)
Edge_KNN5W_Plus$Weight_Reg %>% hist_def(main="W_Reg", save=T)
Edge_KNN5W_Plus$Weight_Corr %>% hist_def(main="W_Corr", save=T)
Edge_KNN5W_Plus$Weight_OVL %>% hist_def(main="W_OVL", save=T)

Edge_KNN5W_Plus_$Weight_STR9 %>% hist_def(main="W_STR9_Scaled", save=T)
Edge_KNN5W_Plus_$Weight_Reg %>% hist_def(main="W_Reg_Scaled", save=T)



dir = "../../processed_data/net_data/BIOCARTA"

# file = sprintf("%s/RNA_Corr_KNN%s+.csv", dir, c(3, 5, 7))
# write.csv(RNA_Corr_KNN3, file=file[1], row.names=F)
# write.csv(RNA_Corr_KNN5, file=file[2], row.names=F)
# write.csv(RNA_Corr_KNN7, file=file[3], row.names=F)

# file = sprintf("%s/KNN%s_STRING_700+.csv", dir, c(3, 5, 7))
# write.csv(STRING_700_KNN3, file=file[1], row.names=F)
# write.csv(STRING_700_KNN5, file=file[2], row.names=F)
# write.csv(STRING_700_KNN7, file=file[3], row.names=F)

# file = sprintf("%s/KNN%s_STRING_900+.csv", dir, c(3, 5, 7))
# write.csv(STRING_900_KNN3, file=file[1], row.names=F)
# write.csv(STRING_900_KNN5, file=file[2], row.names=F)
# write.csv(STRING_900_KNN7, file=file[3], row.names=F)

# file = sprintf("%s/KNN%s_RegNetwork+.csv", dir, c(3, 5, 7))
# write.csv(RegNetwork_KNN3, file=file[1], row.names=F)
# write.csv(RegNetwork_KNN5, file=file[2], row.names=F)
# write.csv(RegNetwork_KNN7, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_STR9_Reg_Corr+.csv", dir, c(3, 5, 7))
# write.csv(Edge_KNN3, file=file[1], row.names=F)
write.csv(Edge_KNN5W, file=file[2], row.names=F)
# write.csv(Edge_KNN7, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_STR9_Reg_Corr+1.csv", dir, c(3, 5, 7))
# write.csv(Edge_KNN3W_Plus, file=file[1], row.names=F)
write.csv(Edge_KNN5W_Plus, file=file[2], row.names=F)
# write.csv(Edge_KNN7W_Plus, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_STR9_Reg_Corr+2.csv", dir, c(3, 5, 7))
# write.csv(Edge_KNN3W_Plus_, file=file[1], row.names=F)
write.csv(Edge_KNN5W_Plus_, file=file[2], row.names=F)
# write.csv(Edge_KNN7W_Plus_, file=file[3], row.names=F)




supplementary = T
if (supplementary) {
  
  save_for_nc = function(df_list, dir=".", num=1, num_fig=NULL, rowNames=F, suppl=T) {
    
    suppressMessages(library(openxlsx))
    is_list = inherits(df_list, "list")
    if (is_list & is.null(num_fig)) num_fig = letters[1:length(df_list)]
    
    if (!suppl) {
      sheets = sprintf("Fig. %s%s", num, num_fig)
      if (!is_list) sheets = sprintf("Fig. %s", num)
      file = sprintf("%s/SourceData_Fig%s.xlsx", dir, num)
    } else {
      sheets = sprintf("Supplementary Fig. %s%s", num, num_fig)
      if (!is_list) sheets = sprintf("Supplementary Fig. %s", num)
      file = sprintf("%s/SourceData_SupplementaryFig%s.xlsx", dir, num)
    }
    
    write.xlsx(df_list, file=file, sheetName=sheets, rowNames=rowNames)
  }
  
  
  ### [Source Data] Supplementary Fig. 1
  STRING_700_SS_ = STRING_700_SS %>% 
    rename(Separation_Score=Dist, Dist_Pathway1=Dist_P1, Dist_Pathway2=Dist_P2, 
           Dist_Pathway12=Dist_P1_P2, Num_Genes_Pathway1=Num1, Num_Genes_Pathway2=Num2, 
           Num_Genes_Overlap=Num_Overlap, Ratio_Genes_Overlap=Ratio_Overlap)
  
  STRING_900_SS_ = STRING_900_SS %>% 
    rename(Separation_Score=Dist, Dist_Pathway1=Dist_P1, Dist_Pathway2=Dist_P2, 
           Dist_Pathway12=Dist_P1_P2, Num_Genes_Pathway1=Num1, Num_Genes_Pathway2=Num2, 
           Num_Genes_Overlap=Num_Overlap, Ratio_Genes_Overlap=Ratio_Overlap)
  
  RegNetwork_SS_ = RegNetwork_SS %>% 
    rename(Separation_Score=Dist, Dist_Pathway1=Dist_P1, Dist_Pathway2=Dist_P2, 
           Dist_Pathway12=Dist_P1_P2, Num_Genes_Pathway1=Num1, Num_Genes_Pathway2=Num2, 
           Num_Genes_Overlap=Num_Overlap, Ratio_Genes_Overlap=Ratio_Overlap)
  
  KNN5_Sim_ = KNN5_Sim %>% rename(Network_X=Net1, Network_Y=Net2)
  PPA_Info_ = list(STRING_700_SS_, STRING_900_SS_, RegNetwork_SS_, KNN5_Sim_)
  PPA_Info_ %>% save_for_nc(num=1, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 2
  STRING_700_KNN5_ = STRING_700_KNN5 %>% 
    subset(select=c(Pathway1, Pathway2, Dist, Ratio_Overlap)) %>% 
    rename(Separation_Score=Dist, Ratio_Genes_Overlap=Ratio_Overlap)
  
  STRING_900_KNN5_ = STRING_900_KNN5 %>% 
    subset(select=c(Pathway1, Pathway2, Dist, Ratio_Overlap)) %>% 
    rename(Separation_Score=Dist, Ratio_Genes_Overlap=Ratio_Overlap)
  
  RegNetwork_KNN5_ = RegNetwork_KNN5 %>% 
    subset(select=c(Pathway1, Pathway2, Dist, Ratio_Overlap)) %>% 
    rename(Separation_Score=Dist, Ratio_Genes_Overlap=Ratio_Overlap)
  
  KNN5_List_ = list(STRING_700_KNN5_, STRING_900_KNN5_, RegNetwork_KNN5_, RNA_Corr_KNN5)
  KNN5_List_ %>% save_for_nc(num=2, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 32
  Corr_List %>% save_for_nc(num=32, suppl=T)
  
  
  ### [Source Data] Fig. 5
  PCA_List = list(PCA_Raw, PCA_Asinh, PCA_ZNorm, PCA_ZGene, 
                  PCA_ssGSEA, PCA_ssGSEA_NN, PCA_GSVA, PCA_SING, PCA_SING_ST)
  
  col = c("Database", "PC1", "PC2")
  select_ = function(df) {
    df = df[, col] %>% mutate(Cell=rownames(.)) %>% relocate(Cell, .before=everything())
    rownames(df) = NULL
    return(df)
  }
  
  PCA_List = PCA_List %>% lapply(select_)
  PCA_List %>% save_for_nc(num=5, suppl=F)
  
  
  ### [Supplementary Data] Supplementary Data 2
  
}