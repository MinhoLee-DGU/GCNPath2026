#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(ggpubr))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))

source("../functions.R")
source("functions.R")
loadings()

options(dplyr.summarise.inform=F)



##### 2. Get Performance Files [Overall]

is_inf_def = function(x) {
  inf_p = sum(x>0 & is.infinite(x))
  inf_m = sum(x<0 & is.infinite(x))
  if (inf_p!=0) sprintf("# Infinite (+) : %s", inf_p) %>% print
  if (inf_m!=0) sprintf("# Infinite (-) : %s", inf_m) %>% print
}

dir1 = sprintf("IC50_GDSC%s", c("", 1, 2))
dir2 = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")

Dir_List = expand.grid(dir1, dir2)
dir_list = Dir_List %>% apply(1, function(x) paste0(x, collapse="/"))

pattern = "pred_test_([0-9]+).csv"


# BMTMKL [O]
dir = "../BMTMKL/Results"
dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_BMTMKL = read_pred_all(dir_list_, "BMTMKL", pattern)

# RF [O]
dir = "../RF/Results"
dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_RF = read_pred_all(dir_list_, "RF", pattern)

# tCNNS [O]
dir = "../tCNNS/Results"
dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_tCNNS = read_pred_all(dir_list_, "tCNNS", pattern)
# "Prediction Inf found in GDSC x Normal... [n=37180 in Fold 4]"
# "Prediction Inf found in GDSC2 x Normal... [n=19828 in Fold 3]"
# "Prediction Inf found in GDSC x Drug_Blind... [n=36970 in Fold 6]"
# "Prediction Inf found in GDSC2 x Drug_Blind... [n=36899 in Fold 6 & 8]"
# "Prediction Inf found in GDSC x Strict_Blind... [n=15086 in Fold 17]"
# "Prediction Inf found... [n=145963 in Fold 3 & 4 & 6 & 8 & 17]"

# Pred_tCNNS$Perf$RMSE %>% is.na %>% sum   # 6
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC" & Test_Type=="Normal" & Fold==4) %>% pull(Prediction) %>% is_inf_def
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC" & Test_Type=="Drug_Blind" & Fold==6) %>% pull(Prediction) %>% is_inf_def
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC" & Test_Type=="Strict_Blind" & Fold==17) %>% pull(Prediction) %>% is_inf_def
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC2" & Test_Type=="Normal" & Fold==3) %>% pull(Prediction) %>% is_inf_def
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC2" & Test_Type=="Drug_Blind" & Fold==6) %>% pull(Prediction) %>% is_inf_def
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC2" & Test_Type=="Drug_Blind" & Fold==8) %>% pull(Prediction) %>% is_inf_def

# HiDRA [O]
dir = "../HiDRA/Results"
dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_HiDRA = read_pred_all(dir_list_, "HiDRA", pattern)

# PaccMann [O]
dir = "../PaccMann/Results"

dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_PaccMann = read_pred_all(dir_list_, "PaccMann", pattern)

# PaccMann_SANGER [O]
dir = "../PaccMann_SANGER/Results"
dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_PaccMann_SG = read_pred_all(dir_list_, "PaccMann_SG", pattern)

# GraphDRP [O]
dir = "../GraphDRP/Results"
dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_GraphDRP = read_pred_all(dir_list_, "GraphDRP", pattern)

# TGDRP & TGSA [O/O]
dir = "../TGSA/Results"
dir_list_tgdrp = sprintf("%s/%s/TGDRP", dir, dir_list)
dir_list_tgsa = sprintf("%s/%s/TGSA", dir, dir_list)
Pred_TGDRP = read_pred_all(dir_list_tgdrp, "TGDRP", pattern)
Pred_TGSA = read_pred_all(dir_list_tgsa, "TGSA", pattern)

# TGDRP_SANGER & TGSA_SANGER [O/O]
dir = "../TGSA_SANGER/Results"
dir_list_tgdrp = sprintf("%s/%s/TGDRP", dir, dir_list)
dir_list_tgsa = sprintf("%s/%s/TGSA", dir, dir_list)
Pred_TGDRP_SG = read_pred_all(dir_list_tgdrp, "TGDRP_SG", pattern)
Pred_TGSA_SG = read_pred_all(dir_list_tgsa, "TGSA_SG", pattern)

# DRPreter & SA [O/O]
dir = "../DRPreter/Results"
dir_list_drp = sprintf("%s/%s/DRPreter", dir, dir_list)
dir_list_drpsa = sprintf("%s/%s/DRPreter_SA", dir, dir_list)
Pred_DRPreter = read_pred_all(dir_list_drp, "DRPreter", pattern)
Pred_DRPreter_SA = read_pred_all(dir_list_drpsa, "DRPreter_SA", pattern)

# DRPreter & SA [O/△]
dir = "../DRPreter_SANGER/Results"
dir_list_drp = sprintf("%s/%s/DRPreter", dir, dir_list)
dir_list_drpsa = sprintf("%s/%s/DRPreter_SA", dir, dir_list)
Pred_DRPreter_SG = read_pred_all(dir_list_drp, "DRPreter_SG", pattern)
Pred_DRPreter_SA_SG = read_pred_all(dir_list_drpsa, "DRPreter_SA_SG", pattern)

# GCNPath [O]
dir = "../GCNPath/results"
dir_list_ = sprintf("%s/%s/RGCN", dir, dir_list)
Pred_GCNPath = read_pred_all(dir_list_, "GCNPath", pattern)

# cf. Examine cell lines
cells_sanger = Pred_GCNPath$Perf_Cell$Cell %>% unique      # 972
cells_drpreter = Pred_DRPreter$Perf_Cell$Cell %>% unique   # 696
cells_tgsa = Pred_TGSA$Perf_Cell$Cell %>% unique           # 700
cells_paccmann = Pred_PaccMann$Perf_Cell$Cell %>% unique   # 395
cells_graphdrp = Pred_GraphDRP$Perf_Cell$Cell %>% unique   # 969
cells_hidra = Pred_HiDRA$Perf_Cell$Cell %>% unique         # 947
cells_rf = Pred_RF$Perf_Cell$Cell %>% unique               # 972
cells_bmtmkl = Pred_BMTMKL$Perf_Cell$Cell %>% unique       # 808

all(cells_drpreter %in% cells_sanger)   # T
all(cells_tgsa %in% cells_sanger)       # T
all(cells_paccmann %in% cells_sanger)   # F [394/395]
all(cells_graphdrp %in% cells_sanger)   # F [964/969]
all(cells_hidra %in% cells_sanger)      # F [944/947]
all(cells_rf %in% cells_sanger)         # T



##### 3-1. Compare Performances [Overall]

model_names = c("BMTMKL", "RF", "tCNNS", "HiDRA", "PaccMann", "PaccMann_SG",
                "GraphDRP", "TGDRP", "TGSA", "TGDRP_SG", "TGSA_SG",
                "DRPreter", "DRPreter_SA", "DRPreter_SG", "DRPreter_SA_SG", "GCNPath")

Perf_List = list()
pred_names = sprintf("Pred_%s", model_names)

for (i in 1:length(model_names)) Perf_List[[i]] = get(pred_names[i])$Perf
Perf_List = Reduce(rbind_perf, Perf_List)   # 2475 x 11

dataset = c("GDSC", "GDSC1", "GDSC2")
test_type = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")

Perf_List = Perf_List %>% 
  mutate(Model = Model %>% factor(levels=model_names), 
         Dataset = Dataset %>% factor(levels=dataset), 
         Test_Type = Test_Type %>% factor(levels=test_type))

Test = expand.grid(dataset, test_type) %>% arrange(Var2)
test = sprintf("%s x %s", Test$Var1, Test$Var2)

# RMSE Performances [Mean & SD]
Perf_List$Test = paste(Perf_List$Dataset, Perf_List$Test_Type, sep=" x ")
Perf_List$Test = Perf_List$Test %>% factor(levels=test)

Perf_RMSE_Avg = Perf_List %>% perf_avg_sd(stat="RMSE")
Perf_PCC_Avg = Perf_List %>% perf_avg_sd(stat="PCC")
Perf_SCC_Avg = Perf_List %>% perf_avg_sd(stat="SCC")

# The number of IC50, Cell, Drug
IC50_Number = Perf_List %>% subset(Test_Type=="Normal") %>% 
  group_by(Model, Dataset) %>% summarise(n=sum(N_Test)) %>% 
  as.data.frame %>% acast(Model~Dataset, value.var="n") %>% as.data.frame

Perf_Cell = list()
Perf_Drug = list()

for (i in 1:length(pred_names)) Perf_Cell[[i]] = get(pred_names[i])$Perf_Cell
for (i in 1:length(pred_names)) Perf_Drug[[i]] = get(pred_names[i])$Perf_Drug

Perf_Cell = Reduce(rbind_perf, Perf_Cell)   # 142028 x 8
Perf_Drug = Reduce(rbind_perf, Perf_Drug)   # 55384 x 8

# Number of cells, drugs
ulen = function(x) x %>% unique %>% length
Perf_Cell$Model = Perf_Cell$Model %>% factor(levels=model_names)
Perf_Drug$Model = Perf_Drug$Model %>% factor(levels=model_names)

Cell_Number = Perf_Cell %>% group_by(Model, Dataset) %>% summarise(n=ulen(Cell)) %>% 
  as.data.frame %>% acast(Model~Dataset, value.var="n") %>% as.data.frame
Drug_Number = Perf_Drug %>% group_by(Model, Dataset) %>% summarise(n=ulen(Drug)) %>% 
  as.data.frame %>% acast(Model~Dataset, value.var="n") %>% as.data.frame

colnames(Cell_Number) = sprintf("Cell_%s", colnames(Cell_Number))
colnames(Drug_Number) = sprintf("Drug_%s", colnames(Drug_Number))
colnames(IC50_Number) = sprintf("IC50_%s", colnames(IC50_Number))

identical(rownames(Cell_Number), rownames(Drug_Number))   # T
identical(rownames(Cell_Number), rownames(IC50_Number))   # T
Stat_Number = Reduce(cbind, list(Cell_Number, Drug_Number, IC50_Number))
Stat_Number = Stat_Number[, c(1, 4, 7, 2, 5, 8, 3, 6, 9)]

file = "Performance_Summary.xlsx"
sheet = c("RMSE_Test", "PCC_Test", "SCC_Test", "Number")
df_list = list(Perf_RMSE_Avg, Perf_PCC_Avg, Perf_SCC_Avg, Stat_Number)
write.xlsx(df_list, file=file, rowNames=T, sheetName=sheet)


# Prediction
save_pred = T
if (save_pred) {
  Pred_List = list()
  for (i in 1:length(model_names)) {
    Pred_List[[i]] = get(pred_names[i])$Pred[Dataset=="GDSC"]
    setkey(Pred_List[[i]], Test_Type, Cell, Drug)
  }
  names(Pred_List) = model_names
  
  dir = "../GCNPath/data/ic50_data"
  file = sprintf("%s/IC50_GDSC.txt", dir)
  IC50_GDSC_ = fread(file, sep="\t")
  
  setkey(IC50_GDSC_, Cell, Drug)
  IC50_GDSC_ = IC50_GDSC_[, .(Cell, Drug, LN_IC50)]
  
  pred_to_wide = function(models, test_type, IC50_GDSC=NULL) {
    long_dt <- rbindlist(
      lapply(models, function(model) {
        dt <- Pred_List[[model]][.(test_type), .(Cell, Drug, Prediction)]
        dt[, Model := model]
        dt
      }), use.names = T
    )
    
    setkey(long_dt, Cell, Drug)
    long_dt = dcast(
      long_dt,
      Cell + Drug ~ Model,
      value.var = "Prediction"
    )
    
    setcolorder(long_dt, c("Cell", "Drug", models))
    if (!is.null(IC50_GDSC)) {
      long_dt = IC50_GDSC[long_dt]
    }
    
    return(long_dt)
  }
  
  model_names_ = model_names %>% setdiff("BMTMKL")
  Pred_U = pred_to_wide(model_names, test_type[1], IC50_GDSC_)
  Pred_C = pred_to_wide(model_names, test_type[2], IC50_GDSC_)
  Pred_D = pred_to_wide(model_names_, test_type[3], IC50_GDSC_)
  Pred_S = pred_to_wide(model_names_, test_type[4], IC50_GDSC_)
  
  # Pred_List = Reduce(rbind_perf, Pred_List)   # 40848552 x 8
  # Pred_List = Pred_List %>% rename(Train_Fold=Fold)
  # 
  # Pred_List_ = rbind(Pred_tCNNS$Pred, Pred_GraphDRP$Pred)
  # Pred_List_ = Pred_List_ %>% rename(Train_Fold=Fold) %>% 
  #   subset(select=-c(Cell_BROAD, Cell_COSMIC))
  # 
  # file = "Prediction [GDSC].csv"
  # fwrite(Pred_List, file=file, row.names=F)
  # 
  # file = "Prediction [GDSC (tCNNS, GraphDRP)].csv"
  # fwrite(Pred_List_, file=file, row.names=F)
  # 
  # rm(Pred_List, Pred_List_)
  # gc()
}


# Scatter plots of prediction-actual ln(IC50) values
# Be careful that it takes too much times to draw all figures...

draw_pred = T
if (draw_pred) {
  Info_N = data.frame()
  Info_C = data.frame()
  Info_D = data.frame()
  Info_S = data.frame()
  
  dir_n = mkdir("Prediction/Normal")
  dir_c = mkdir("Prediction/Cell_Blind")
  dir_d = mkdir("Prediction/Drug_Blind")
  dir_s = mkdir("Prediction/Strict_Blind")
  
  for (i in 1:length(pred_names)) {
    Pred_Temp = get(pred_names[i])$Pred
    model = gsub("Pred_", "", pred_names[i])
    
    info_n = Pred_Temp %>% 
      subset(Test_Type=="Normal" & Dataset=="GDSC") %>% 
      plot_pred(dir=dir_n, model=model, test_type="Normal")
    info_c = Pred_Temp %>% 
      subset(Test_Type=="Cell_Blind" & Dataset=="GDSC") %>% 
      plot_pred(dir=dir_c, model=model, test_type="Cell-Blind")
    info_d = Pred_Temp %>% 
      subset(Test_Type=="Drug_Blind" & Dataset=="GDSC") %>% 
      plot_pred(dir=dir_d, model=model, test_type="Drug-Blind")
    info_s = Pred_Temp %>% 
      subset(Test_Type=="Strict_Blind" & Dataset=="GDSC") %>% 
      plot_pred(dir=dir_s, model=model, test_type="Strict-Blind")
    
    Info_N = Info_N %>% rbind(c(model, info_n))
    Info_C = Info_C %>% rbind(c(model, info_c))
    Info_D = Info_D %>% rbind(c(model, info_d))
    Info_S = Info_S %>% rbind(c(model, info_s))
  }
  
  col = c("Model", "N_Test", "RMSE", "PCC")
  Info_N = Info_N %>% setNames(col)
  Info_C = Info_C %>% setNames(col)
  Info_D = Info_D %>% setNames(col)
  Info_S = Info_S %>% setNames(col)
}


# Statistic Analysis per Dataset & Test-Type
# Non-parametric U-Test [Mann-Whitney]
Perf_UTest = Perf_List %>% wilcox_test(fdr_adjust=T)

# Confirmation of FDR adjust p.val
utest_is_fdr = c()
for (dataset_ in dataset) {
  for (test_type_ in test_type) {
    Temp = Perf_UTest %>% subset(Dataset==dataset_ & Test_Type==test_type_)
    utest_is_fdr_ = identical(Temp$FDR_RMSE, p.adjust(Temp$Pval_RMSE, "fdr"))
    utest_is_fdr = utest_is_fdr %>% c(utest_is_fdr_)
  }
}

all(utest_is_fdr)   # T
test_type_ = gsub("_", "-", test_type)
Perf_UTest_N = Perf_UTest %>% subset(Dataset=="GDSC" & Test_Type==test_type[1])
Perf_UTest_C = Perf_UTest %>% subset(Dataset=="GDSC" & Test_Type==test_type[2])
Perf_UTest_D = Perf_UTest %>% subset(Dataset=="GDSC" & Test_Type==test_type[3])
Perf_UTest_S = Perf_UTest %>% subset(Dataset=="GDSC" & Test_Type==test_type[4])

dir = mkdir("Performance [Wilcox Test, Grid]")
main1 = sprintf("%s/%s Wilcox of Models [GDSC, %s]", dir, metrics, test_type_[1])
main2 = sprintf("%s/%s Wilcox of Models [GDSC, %s]", dir, metrics, test_type_[2])
main3 = sprintf("%s/%s Wilcox of Models [GDSC, %s]", dir, metrics, test_type_[3])
main4 = sprintf("%s/%s Wilcox of Models [GDSC, %s]", dir, metrics, test_type_[4])

model_names_ = setdiff(model_names, "BMTMKL")
Perf_UTest_N %>% wilcox_model_grid(models=model_names, lvl_model=model_names, main=main1, save=T)
Perf_UTest_C %>% wilcox_model_grid(models=model_names, lvl_model=model_names, main=main2, save=T)
Perf_UTest_D %>% wilcox_model_grid(models=model_names_, lvl_model=model_names_, main=main3, save=T)
Perf_UTest_S %>% wilcox_model_grid(models=model_names_, lvl_model=model_names_, main=main4, save=T)

dir = mkdir("Performance [Overall, Boxplot]")

Perf_List %>% plot_perf("RMSE", "Normal", dir=dir, save=T)
Perf_List %>% plot_perf("RMSE", "Cell_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("RMSE", "Drug_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("RMSE", "Strict_Blind", dir=dir, save=T)

Perf_List %>% plot_perf("PCC", "Normal", dir=dir, save=T)
Perf_List %>% plot_perf("PCC", "Cell_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("PCC", "Drug_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("PCC", "Strict_Blind", dir=dir, save=T)

Perf_List %>% plot_perf("SCC", "Normal", dir=dir, save=T)
Perf_List %>% plot_perf("SCC", "Cell_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("SCC", "Drug_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("SCC", "Strict_Blind", dir=dir, save=T)

dir = mkdir("Performance [Overall, Boxplot (Except)]")
models_except = c("tCNNS", "HiDRA", "GraphDRP", "PaccMann", "PaccMann_SG", "RF")

Perf_List %>% 
  subset(!(Model %in% models_except)) %>% 
  mutate(Model=droplevels(Model)) %>%
  plot_perf("RMSE", "Normal", dir=dir, width=30, save=T)

Perf_List %>% 
  subset(!(Model %in% models_except)) %>% 
  mutate(Model=droplevels(Model)) %>%
  plot_perf("PCC", "Normal", dir=dir, width=36, save=T)

Perf_List %>% 
  subset(!(Model %in% models_except)) %>% 
  mutate(Model=droplevels(Model)) %>%
  plot_perf("SCC", "Normal", dir=dir, width=36, save=T)


# Barplot [GDSC1+2 Test]
dir = mkdir("Performance [Overall, Barplot]")
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[1], dir=dir, save=T)
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[2], dir=dir, save=T)
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[3], dir=dir, save=T)
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[4], dir=dir, save=T)

model_except_n = c("tCNNS", "HiDRA", "PaccMann", "PaccMann_SG", "GraphDRP", "RF", "BMTMKL")
model_except_c = c("tCNNS", "GraphDRP")
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[1],
                            model_except=model_except_n, dir=dir, width=20, axis_tx=16.5, save=T)
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[2],
                            model_except=model_except_c, dir=dir, width=20, axis_tx=16.5, save=T)



##### 3-2. Compare Performances [Cell & Drug]

# Cell Annotation from SANGER Cell Passports
file = "Anno_Cells.csv"
Anno_Cells = read.csv(file)

# Drug Annotation from GDSC
file = "Anno_Drugs.csv"
Anno_Drugs = read.csv(file)

idx1 = match(Perf_Cell$Cell, Anno_Cells$SANGER_MODEL_ID)
idx2 = match(Perf_Drug$Drug, Anno_Drugs$Drug_CID)
Perf_Cell = Perf_Cell %>% mutate(Cell_TCGA=Anno_Cells$TCGA_CODE[idx1]) %>% relocate(Cell_TCGA, .after=Cell)
Perf_Drug = Perf_Drug %>% mutate(Drug_Pathway=Anno_Drugs$Target_Pathway[idx2]) %>% relocate(Drug_Pathway, .after=Drug)

Perf_Cell$Cell_TCGA %>% is.na %>% sum   # 0
Perf_Drug$Drug_Pathway %>% is.na %>% sum   # 224
Perf_Cell$Cell_TCGA[is.na(Perf_Cell$Cell_TCGA)] = "UNCLASSIFIED"
Perf_Drug$Drug_Pathway[is.na(Perf_Drug$Drug_Pathway)] = "Unclassified"

Perf_Cell_N = Perf_Cell %>% subset(Dataset=="GDSC" & Test_Type==test_type[1])
Perf_Cell_C = Perf_Cell %>% subset(Dataset=="GDSC" & Test_Type==test_type[2])
Perf_Cell_D = Perf_Cell %>% subset(Dataset=="GDSC" & Test_Type==test_type[3])
Perf_Cell_S = Perf_Cell %>% subset(Dataset=="GDSC" & Test_Type==test_type[4])

Perf_Drug_N = Perf_Drug %>% subset(Dataset=="GDSC" & Test_Type==test_type[1])
Perf_Drug_C = Perf_Drug %>% subset(Dataset=="GDSC" & Test_Type==test_type[2])
Perf_Drug_D = Perf_Drug %>% subset(Dataset=="GDSC" & Test_Type==test_type[3])
Perf_Drug_S = Perf_Drug %>% subset(Dataset=="GDSC" & Test_Type==test_type[4])


dir = mkdir("Performance [Cell and Drug, Grid]")
metrics = c("RMSE", "PCC", "SCC")

test_type_ = gsub("_", "-", test_type)
file1 = sprintf("%s/%s Corr of Models [Cell, %s]", dir, metrics, test_type_[1])
file2 = sprintf("%s/%s Corr of Models [Cell, %s]", dir, metrics, test_type_[2])
file3 = sprintf("%s/%s Corr of Models [Cell, %s]", dir, metrics, test_type_[3])
file4 = sprintf("%s/%s Corr of Models [Cell, %s]", dir, metrics, test_type_[4])
model_names_ = setdiff(model_names, "BMTMKL")

Perf_Cell_Corr_N = Perf_Cell_N %>% 
  compare_model_grid(models=model_names, by="Cell", lvl_model=model_names, main=file1, save=T)
Perf_Cell_Corr_C = Perf_Cell_C %>% 
  compare_model_grid(models=model_names, by="Cell", lvl_model=model_names, main=file2, save=T)
Perf_Cell_Corr_D = Perf_Cell_D %>% 
  compare_model_grid(models=model_names_, by="Cell", lvl_model=model_names_, main=file3, save=T)
Perf_Cell_Corr_S = Perf_Cell_S %>% 
  compare_model_grid(models=model_names_, by="Cell", lvl_model=model_names_, main=file4, save=T)


file1 = sprintf("%s/%s Corr of Models [Drug, %s]", dir, metrics, test_type_[1])
file2 = sprintf("%s/%s Corr of Models [Drug, %s]", dir, metrics, test_type_[2])
file3 = sprintf("%s/%s Corr of Models [Drug, %s]", dir, metrics, test_type_[3])
file4 = sprintf("%s/%s Corr of Models [Drug, %s]", dir, metrics, test_type_[4])

Perf_Drug_Corr_N = Perf_Drug_N %>% 
  compare_model_grid(models=model_names, by="Drug", lvl_model=model_names, main=file1, save=T)
Perf_Drug_Corr_C = Perf_Drug_C %>% 
  compare_model_grid(models=model_names, by="Drug", lvl_model=model_names, main=file2, save=T)
Perf_Drug_Corr_D = Perf_Drug_D %>% 
  compare_model_grid(models=model_names_, by="Drug", lvl_model=model_names_, main=file3, save=T)
Perf_Drug_Corr_S = Perf_Drug_S %>% 
  compare_model_grid(models=model_names_, by="Drug", lvl_model=model_names_, main=file4, save=T)

# Number of cells, drugs
ulen = function(x) x %>% unique %>% length
Perf_Cell$Model = Perf_Cell$Model %>% factor(levels=model_names)
Perf_Drug$Model = Perf_Drug$Model %>% factor(levels=model_names)

Cell_Number = Perf_Cell %>% group_by(Model, Dataset) %>% summarise(n=ulen(Cell)) %>% 
  as.data.frame %>% acast(Model~Dataset, value.var="n") %>% as.data.frame
Drug_Number = Perf_Drug %>% group_by(Model, Dataset) %>% summarise(n=ulen(Drug)) %>% 
  as.data.frame %>% acast(Model~Dataset, value.var="n") %>% as.data.frame

colnames(Cell_Number) = sprintf("Cell_%s", colnames(Cell_Number))
colnames(Drug_Number) = sprintf("Drug_%s", colnames(Drug_Number))
colnames(IC50_Number) = sprintf("IC50_%s", colnames(IC50_Number))

identical(rownames(Cell_Number), rownames(Drug_Number))   # T
identical(rownames(Cell_Number), rownames(IC50_Number))   # T
Stat_Number = Reduce(cbind, list(Cell_Number, Drug_Number, IC50_Number))
Stat_Number = Stat_Number[, c(1, 4, 7, 2, 5, 8, 3, 6, 9)]

Perf_Cell_Corr = list(Perf_Cell_Corr_N, Perf_Cell_Corr_C, Perf_Cell_Corr_D, Perf_Cell_Corr_S)
Perf_Drug_Corr = list(Perf_Drug_Corr_N, Perf_Drug_Corr_C, Perf_Drug_Corr_D, Perf_Drug_Corr_S)

names(Perf_Cell_Corr) = test_type
names(Perf_Drug_Corr) = test_type

file = sprintf("%s/Performance_Cell.csv", dir)
write.csv(Perf_Cell, file=file, row.names=F)

file = sprintf("%s/Performance_Drug.csv", dir)
write.csv(Perf_Drug, file=file, row.names=F)



##### 3-3. Compare Performances [Cell TCGA & Drug MoA]
Pred_GDSC = list()
col = c("Model", "Dataset", "Test_Type", 
        "Fold", "Cell", "Drug", "LN_IC50", "Prediction")

for (i in 1:length(model_names)) {
  Pred_GDSC[[i]] = get(pred_names[i])$Pred[, col, with=F]
}

Pred_GDSC = rbindlist(Pred_GDSC)
idx1 = match(Pred_GDSC$Cell, Anno_Cells$SANGER_MODEL_ID)
idx2 = match(Pred_GDSC$Drug, Anno_Drugs$Drug_CID)
Pred_GDSC[, Cell_TCGA:=Anno_Cells$TCGA_CODE[idx1]]
Pred_GDSC[, Drug_Pathway:=Anno_Drugs$Target_Pathway[idx2]]

Pred_GDSC$Cell_TCGA %>% is.na %>% sum      # 0
Pred_GDSC$Drug_Pathway %>% is.na %>% sum   # 0
Pred_GDSC$Drug_Pathway[is.na(Pred_GDSC$Drug_Pathway)] = "Unclassified"

use = "pairwise.complete.obs"
Perf_Cell_TCGA = Pred_GDSC[!is.infinite(Prediction), 
                           .(RMSE = RMSE(Prediction, LN_IC50, na.rm=T), 
                             PCC = cor(Prediction, LN_IC50, use=use), 
                             SCC = cor(Prediction, LN_IC50, use=use, method="spearman")), 
                           by=.(Model, Dataset, Test_Type, Cell_TCGA)]

Perf_Drug_Pathway = Pred_GDSC[!is.infinite(Prediction), 
                              .(RMSE = RMSE(Prediction, LN_IC50, na.rm=T), 
                                PCC = cor(Prediction, LN_IC50, use=use), 
                                SCC = cor(Prediction, LN_IC50, use=use, method="spearman")), 
                              by=.(Model, Dataset, Test_Type, Drug_Pathway)]

Perf_Cell_TCGA = Perf_Cell_TCGA %>% as.data.frame
Perf_Drug_Pathway = Perf_Drug_Pathway %>% as.data.frame

dataset_ = "GDSC2"

### Cell TCGA & Drug MoA [Barplot]
dir = mkdir(sprintf("Performance by Cell TCGA and Drug MoA [Barplot, %s]", dataset_))
main = sprintf("%s/Performance by Cell TCGA [%s]", dir, test_type)
Perf_Cell_TCGA %>% plot_perf_bar(score="RMSE", dataset=dataset_, test_type=test_type[1], main=main[1], save=T)
Perf_Cell_TCGA %>% plot_perf_bar(score="RMSE", dataset=dataset_, test_type=test_type[2], main=main[2], save=T)
Perf_Cell_TCGA %>% plot_perf_bar(score="RMSE", dataset=dataset_, test_type=test_type[4], main=main[4], save=T)

main = sprintf("%s/Performance by Drug MoA [%s]", dir, test_type)
Perf_Drug_Pathway %>% plot_perf_bar(score="RMSE", dataset=dataset_, test_type=test_type[1], main=main[1], save=T)
Perf_Drug_Pathway %>% plot_perf_bar(score="RMSE", dataset=dataset_, test_type=test_type[3], main=main[3], save=T)
Perf_Drug_Pathway %>% plot_perf_bar(score="RMSE", dataset=dataset_, test_type=test_type[4], main=main[4], save=T)

### Cell TCGA & Drug MoA [Heatmap]
dir = mkdir("Performance by Cell TCGA and Drug MoA [Heatmap, GDSC1+2]")
main = sprintf("%s/Performance by Cell TCGA [%s]", dir, test_type)

Perf_Cell_TCGA %>% 
  subset(Dataset=="GDSC" & Test_Type==test_type[1]) %>% 
  mutate(Model=factor(Model, levels=model_names)) %>% 
  reshape2::acast(Cell_TCGA~Model, value.var="RMSE") %>% as.data.frame %>% 
  heatmap_save(main=main[1], add=add, width=24, height=40, text_num=12, 
               show_row=T, show_col=T, show_num=T, center_zero=T, 
               clust_row=F, clust_col=F, scale_row=T, save=T, save_svg=T)

Perf_Cell_TCGA %>% 
  subset(Dataset=="GDSC" & Test_Type==test_type[2]) %>% 
  mutate(Model=factor(Model, levels=model_names)) %>% 
  reshape2::acast(Cell_TCGA~Model, value.var="RMSE") %>% as.data.frame %>% 
  heatmap_save(main=main[2], add=add, width=24, height=40, text_num=12, 
               show_row=T, show_col=T, show_num=T, center_zero=T, 
               clust_row=F, clust_col=F, scale_row=T, save=T, save_svg=T)

Perf_Cell_TCGA %>% 
  subset(Dataset=="GDSC" & Test_Type==test_type[4]) %>% 
  mutate(Model=factor(Model, levels=model_names)) %>% 
  reshape2::acast(Cell_TCGA~Model, value.var="RMSE") %>% as.data.frame %>% 
  heatmap_save(main=main[4], add=add, width=24, height=40, text_num=12, 
               show_row=T, show_col=T, show_num=T, center_zero=T, 
               clust_row=F, clust_col=F, scale_row=T, save=T, save_svg=T)

main = sprintf("%s/Performance by Drug Pathway [%s]", dir, test_type)

Perf_Drug_Pathway %>% 
  subset(Dataset=="GDSC" & Test_Type==test_type[1]) %>% 
  mutate(Model=factor(Model, levels=model_names)) %>% 
  reshape2::acast(Drug_Pathway~Model, value.var="RMSE") %>% as.data.frame %>% 
  heatmap_save(main=main[1], add=add, width=27, height=30, text_num=12, 
               show_row=T, show_col=T, show_num=T, center_zero=T, 
               clust_row=F, clust_col=F, scale_row=T, save=T, save_svg=T)

Perf_Drug_Pathway %>% 
  subset(Dataset=="GDSC" & Test_Type==test_type[3]) %>% 
  mutate(Model=factor(Model, levels=model_names)) %>% 
  reshape2::acast(Drug_Pathway~Model, value.var="RMSE") %>% as.data.frame %>% 
  heatmap_save(main=main[3], add=add, width=27, height=30, text_num=12, 
               show_row=T, show_col=T, show_num=T, center_zero=T, 
               clust_row=F, clust_col=F, scale_row=T, save=T, save_svg=T)

Perf_Drug_Pathway %>% 
  subset(Dataset=="GDSC" & Test_Type==test_type[4]) %>% 
  mutate(Model=factor(Model, levels=model_names)) %>% 
  reshape2::acast(Drug_Pathway~Model, value.var="RMSE") %>% as.data.frame %>% 
  heatmap_save(main=main[4], add=add, width=27, height=30, text_num=12, 
               show_row=T, show_col=T, show_num=T, center_zero=T, 
               clust_row=F, clust_col=F, scale_row=T, save=T, save_svg=T)


add1 = list(scale_fill_manual(values=rep("#66b3ed", 38)))
add2 = list(scale_fill_manual(values=rep("#66b3ed", 24)))
dir = mkdir("Performance by Cell TCGA and Drug MoA [GCNPath]")

ylab = bquote(RMSE[C])
main = sprintf("%s/Performance by Cell TCGA in GCNPath [GDSC1+2, %s]", dir, test_type)
Perf_Cell %>% 
  subset(Model=="GCNPath" & Dataset=="GDSC" & Test_Type==test_type[1]) %>%
  boxplot_def(Cell_TCGA, RMSE, main=main[1], legend=F, reorder=T, axis_tl=24, add=add1, alpha=0.9,
              ylab=ylab, angle=36, width=32, hjust=1, vjust=1, legend_tx=10.8, save=T)
Perf_Cell %>% 
  subset(Model=="GCNPath" & Dataset=="GDSC" & Test_Type==test_type[2]) %>%
  boxplot_def(Cell_TCGA, RMSE, main=main[2], legend=F, reorder=T, axis_tl=24, add=add1, alpha=0.9,
              ylab=ylab, angle=36, width=32, hjust=1, vjust=1, legend_tx=10.8, save=T)
Perf_Cell %>% 
  subset(Model=="GCNPath" & Dataset=="GDSC" & Test_Type==test_type[4]) %>%
  boxplot_def(Cell_TCGA, RMSE, main=main[4], legend=F, reorder=T, axis_tl=24, add=add1, alpha=0.9,
              ylab=ylab, angle=36, width=32, hjust=1, vjust=1, legend_tx=10.8, save=T)

ylab = bquote(RMSE[D])
main = sprintf("%s/Performance by Drug MoA in GCNPath [GDSC1+2, %s]", dir, test_type)
Perf_Drug %>% 
  subset(Model=="GCNPath" & Dataset=="GDSC" & Test_Type==test_type[1]) %>%
  boxplot_def(Drug_Pathway, RMSE, main=main[1], legend=F, reorder=T, axis_tl=24, add=add2, alpha=0.9,
              ylab=ylab, angle=36, width=32, hjust=1, vjust=1, legend_tx=9, save=T)
Perf_Drug %>% 
  subset(Model=="GCNPath" & Dataset=="GDSC" & Test_Type==test_type[3]) %>%
  boxplot_def(Drug_Pathway, RMSE, main=main[3], legend=F, reorder=T, axis_tl=24, add=add2, alpha=0.9,
              ylab=ylab, angle=36, width=32, hjust=1, vjust=1, legend_tx=9, save=T)
Perf_Drug %>% 
  subset(Model=="GCNPath" & Dataset=="GDSC" & Test_Type==test_type[4]) %>%
  boxplot_def(Drug_Pathway, RMSE, main=main[4], legend=F, reorder=T, axis_tl=24, add=add2, alpha=0.9,
              ylab=ylab, angle=36, width=32, hjust=1, vjust=1, legend_tx=9, save=T)




# [GCNPath, GDSC1 vs GDSC1+2]
dataset = c("GDSC", "GDSC1", "GDSC2")
dir = mkdir("Performance by Cell TCGA and Drug MoA [GCNPath, GDSC1 vs GDSC1+2]")
main1 = sprintf("%s/Performance by Cell TCGA in GCNPath [%s]", dir, test_type[c(1, 2, 4)])
main2 = sprintf("%s/Performance by Drug MoA in GCNPath [%s]", dir, test_type[c(1, 3, 4)])

Perf_GCNPath_GDSC1C = Perf_Cell_TCGA %>% 
  subset(Model=="GCNPath" & Dataset %in% dataset[c(1, 2)] & Test_Type %in% test_type[c(1, 2, 4)]) %>%
  reshape2::dcast(Cell_TCGA~Dataset+Test_Type, value.var="RMSE") %>% 
  relocate(GDSC1_Normal, GDSC1_Cell_Blind, GDSC1_Strict_Blind, 
           GDSC_Normal, GDSC_Cell_Blind, GDSC_Strict_Blind, .after=Cell_TCGA)

Perf_GCNPath_GDSC1D = Perf_Drug_Pathway %>% 
  subset(Model=="GCNPath" & Dataset %in% dataset[c(1, 2)] & Test_Type %in% test_type[c(1, 3, 4)]) %>%
  reshape2::dcast(Drug_Pathway~Dataset+Test_Type, value.var="RMSE") %>% 
  relocate(GDSC1_Normal, GDSC1_Drug_Blind, GDSC1_Strict_Blind, 
           GDSC_Normal, GDSC_Drug_Blind, GDSC_Strict_Blind, .after=Drug_Pathway)

add_c = list(text_repel_def(NULL, Cell_TCGA))
add_d = list(text_repel_def(NULL, Drug_Pathway))
labs = sprintf("RMSE [%s]", dataset[c(1, 2)])

Perf_GCNPath_GDSC1C %>% 
  plot_def(GDSC1_Normal, GDSC_Normal, 
           main=main1[1], xlab=labs[1], ylab=labs[2], add=add_c, 
           width=10, height=10, xy_line=T, unify_lims=T, save=T)

Perf_GCNPath_GDSC1C %>% 
  plot_def(GDSC1_Cell_Blind, GDSC_Cell_Blind, 
           main=main1[2], xlab=labs[1], ylab=labs[2], add=add_c, 
           width=10, height=10, xy_line=T, unify_lims=T, save=T)

Perf_GCNPath_GDSC1C %>% 
  plot_def(GDSC1_Strict_Blind, GDSC_Strict_Blind, 
           main=main1[3], xlab=labs[1], ylab=labs[2], add=add_c, 
           width=10, height=10, xy_line=T, unify_lims=T, save=T)

Perf_GCNPath_GDSC1D %>% 
  plot_def(GDSC1_Normal, GDSC_Normal, 
           main=main2[1], xlab=labs[1], ylab=labs[2], add=add_d, 
           width=10, height=10, xy_line=T, unify_lims=T, save=T)

Perf_GCNPath_GDSC1D %>% 
  plot_def(GDSC1_Drug_Blind, GDSC_Drug_Blind, 
           main=main2[2], xlab=labs[1], ylab=labs[2], add=add_d,  
           width=10, height=10, xy_line=T, unify_lims=T, save=T)

Perf_GCNPath_GDSC1D %>% 
  plot_def(GDSC1_Strict_Blind, GDSC_Strict_Blind, 
           main=main2[3], xlab=labs[1], ylab=labs[2], add=add_d, 
           width=10, height=10, xy_line=T, unify_lims=T, save=T)


# [GCNPath, GDSC2 vs GDSC1+2]
labs = sprintf("RMSE [%s]", dataset[c(1, 3)])
dir = mkdir("Performance by Cell TCGA and Drug MoA [GCNPath, GDSC2 vs GDSC1+2]")
main1 = sprintf("%s/Performance by Cell TCGA in GCNPath [%s]", dir, test_type[c(1, 2, 4)])
main2 = sprintf("%s/Performance by Drug MoA in GCNPath [%s]", dir, test_type[c(1, 3, 4)])

Perf_GCNPath_GDSC2C = Perf_Cell_TCGA %>% 
  subset(Model=="GCNPath" & Dataset %in% dataset[c(1, 3)] & Test_Type %in% test_type[c(1, 2, 4)]) %>%
  reshape2::dcast(Cell_TCGA~Dataset+Test_Type, value.var="RMSE") %>% 
  relocate(GDSC2_Normal, GDSC2_Cell_Blind, GDSC2_Strict_Blind, 
           GDSC_Normal, GDSC_Cell_Blind, GDSC_Strict_Blind, .after=Cell_TCGA)

Perf_GCNPath_GDSC2D = Perf_Drug_Pathway %>% 
  subset(Model=="GCNPath" & Dataset %in% dataset[c(1, 3)] & Test_Type %in% test_type[c(1, 3, 4)]) %>%
  reshape2::dcast(Drug_Pathway~Dataset+Test_Type, value.var="RMSE") %>% 
  relocate(GDSC2_Normal, GDSC2_Drug_Blind, GDSC2_Strict_Blind, 
           GDSC_Normal, GDSC_Drug_Blind, GDSC_Strict_Blind, .after=Drug_Pathway)

Perf_GCNPath_GDSC2C %>% 
  plot_def(GDSC2_Normal, GDSC_Normal, 
           main=main1[1], xlab=labs[1], ylab=labs[2], add=add_c,
           width=10, height=10, xy_line=T, unify_lims=T, save=T)

Perf_GCNPath_GDSC2C %>% 
  plot_def(GDSC2_Cell_Blind, GDSC_Cell_Blind, 
           main=main1[2], xlab=labs[1], ylab=labs[2], add=add_c,
           width=10, height=10, xy_line=T, unify_lims=T, save=T)

Perf_GCNPath_GDSC2C %>% 
  plot_def(GDSC2_Strict_Blind, GDSC_Strict_Blind, 
           main=main1[3], xlab=labs[1], ylab=labs[2], add=add_c,
           width=10, height=10, xy_line=T, unify_lims=T, save=T)

Perf_GCNPath_GDSC2D %>% 
  plot_def(GDSC2_Normal, GDSC_Normal, 
           main=main2[1], xlab=labs[1], ylab=labs[2], add=add_d,
           width=10, height=10, xy_line=T, unify_lims=T, save=T)

Perf_GCNPath_GDSC2D %>% 
  plot_def(GDSC2_Drug_Blind, GDSC_Drug_Blind,
           main=main2[2], xlab=labs[1], ylab=labs[2], add=add_d,
           width=10, height=10, xy_line=T, unify_lims=T, save=T)

Perf_GCNPath_GDSC2D %>% 
  plot_def(GDSC2_Strict_Blind, GDSC_Strict_Blind, 
           main=main2[3], xlab=labs[1], ylab=labs[2], add=add_d,
           width=10, height=10, xy_line=T, unify_lims=T, save=T)




##### 4. Performance Examination [Detailed]
### Which factors affect the RMSE (based on GCNPath)?

bimodality <- function(x) {
  x <- na.omit(x)
  n <- length(x)
  g <- e1071::skewness(x, type = 1)
  k <- e1071::kurtosis(x, type = 1) + 3
  bc <- (g^2 + 1) / (k + (3 * (n - 1)^2) / ((n - 2) * (n - 3)))
  return(bc)
}

dcoverage <- function(x, global_cdf) {
  q10 <- quantile(x, 0.10, na.rm = TRUE)
  q90 <- quantile(x, 0.90, na.rm = TRUE)
  return(global_cdf(q90) - global_cdf(q10))
}

# IC50 CCLE
# GDSC_Last/processed_data/ic50_data/IC50_CCLE.csv
IC50_CCLE = read.csv("IC50_CCLE.csv")
IC50_CCLE = IC50_CCLE %>% rename(LN_IC50_CCLE=LN_IC50) %>% subset(!Capping)   # 5171
# IC50_CCLE, IC50_CCLE_FT

Stat_Cell_List = list()
Stat_Drug_List = list()
Stat_CCLE_List = list()
Stat_Cell_List_Group = list()
Stat_Drug_List_Group = list()

add = list(stat_smooth(method="lm"))
test_type = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")

for (ttype in test_type) {
  sprintf("### Performance Analysis in %s", ttype) %>% print
  dir = mkdir(sprintf("Analysis of Error [%s]", ttype))
  
  Pred_GCNPath_ = Pred_GCNPath$Pred %>% as.data.frame %>%
    subset(Test_Type==ttype & Dataset=="GDSC")
  
  # Anno_Cells, Anno_Drugs
  idx1 = match(Pred_GCNPath_$Cell, Anno_Cells$SANGER_MODEL_ID)
  idx2 = match(Pred_GCNPath_$Drug, Anno_Drugs$Drug_CID)
  
  Pred_GCNPath_ = Pred_GCNPath_ %>%
    mutate(Cell_TCGA=Anno_Cells$TCGA_CODE[idx1],
           Drug_Pathway=Anno_Drugs$Target_Pathway[idx2])
  
  Pred_GCNPath_$Cell_TCGA %>% is.na %>% sum      # 0
  Pred_GCNPath_$Drug_Pathway %>% is.na %>% sum   # 1426
  Pred_GCNPath_$Drug_Pathway[is.na(Pred_GCNPath_$Drug_Pathway)] = "Unclassified"
  
  global_cdf = ecdf(Pred_GCNPath_$LN_IC50)
  col = c("Cell", "Drug", "LN_IC50", "Prediction", "Cell_TCGA", "Drug_Pathway")
  
  IC50_Stat_Cell = Pred_GCNPath_[, col] %>% group_by(Cell, Cell_TCGA) %>%
    summarise(Num=n(),
              SCC=cor(LN_IC50, Prediction, method="spearman"),
              IC50_Mean=mean(LN_IC50), IC50_SD=sd(LN_IC50),
              IC50_Bimodality=bimodality(LN_IC50),
              IC50_Coverage=dcoverage(LN_IC50, global_cdf)) %>% as.data.frame   # 972
  
  IC50_Stat_Drug = Pred_GCNPath_[, col] %>% group_by(Drug, Drug_Pathway) %>%
    summarise(Num=n(),
              SCC=cor(LN_IC50, Prediction, method="spearman"),
              IC50_Mean=mean(LN_IC50), IC50_SD=sd(LN_IC50),
              IC50_Bimodality=bimodality(LN_IC50),
              IC50_Coverage=dcoverage(LN_IC50, global_cdf)) %>% as.data.frame   # 432
  
  IC50_Stat_Cell_Group = Pred_GCNPath_[, col] %>% group_by(Cell_TCGA) %>%
    summarise(Num=n(), SCC=cor(LN_IC50, Prediction, method="spearman")) %>% as.data.frame   # 38
  
  IC50_Stat_Drug_Group = Pred_GCNPath_[, col] %>% group_by(Drug_Pathway) %>%
    summarise(Num=n(), SCC=cor(LN_IC50, Prediction, method="spearman")) %>% as.data.frame   # 24
  
  Stat_Cell_List[[ttype]] = IC50_Stat_Cell
  Stat_Drug_List[[ttype]] = IC50_Stat_Drug
  Stat_Cell_List_Group[[ttype]] = IC50_Stat_Cell_Group
  Stat_Drug_List_Group[[ttype]] = IC50_Stat_Drug_Group
  
  
  ### 1. IC50 Variation
  print("# 1. IC50 Variation")
  xlab = bquote(SD~ln(IC[50]))
  main = sprintf("%s/IC50_SD & SCC Relation [%s]", dir, by)
  
  IC50_Stat_Drug %>% with(cor(IC50_SD, SCC)) %>% round(3) %>% print   # PCC=0.776
  IC50_Stat_Drug %>% plot_def(IC50_SD, SCC, main=main[2], xlab=xlab, ylab=ylab_d, add=add,
                              size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)
  
  IC50_Stat_Cell %>% with(cor(IC50_SD, SCC)) %>% round(3) %>% print   # PCC=-0.003
  IC50_Stat_Cell %>% plot_def(IC50_SD, SCC, main=main[1], xlab=xlab, ylab=ylab_c, add=add,
                              size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)
  
  
  ### 2. IC50 Bimodality
  print("# 2. IC50 Bimodality")
  xlab = bquote(Bimodality~ln(IC[50]))
  main = sprintf("%s/IC50_Bimodality & SCC Relation [%s]", dir, by)
  
  IC50_Stat_Drug %>% with(cor(IC50_Bimodality, SCC)) %>% round(3) %>% print   # PCC=0.776
  IC50_Stat_Drug %>% plot_def(IC50_Bimodality, SCC, main=main[2], xlab=xlab, ylab=ylab_d, add=add,
                              size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)
  
  IC50_Stat_Cell %>% with(cor(IC50_Bimodality, SCC)) %>% round(3) %>% print   # PCC=-0.003
  IC50_Stat_Cell %>% plot_def(IC50_Bimodality, SCC, main=main[1], xlab=xlab, ylab=ylab_c, add=add,
                              size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)
  
  
  ### 3. IC50 Coverage
  print("# 3. IC50 Coverage")
  xlab = bquote(Coverage~ln(IC[50]))
  main = sprintf("%s/IC50_Coverage & SCC Relation [%s]", dir, by)
  
  IC50_Stat_Drug %>% with(cor(IC50_Coverage, SCC)) %>% round(3) %>% print   # PCC=0.776
  IC50_Stat_Drug %>% plot_def(IC50_Coverage, SCC, main=main[2], xlab=xlab, ylab=ylab_d, add=add,
                              size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)
  
  IC50_Stat_Cell %>% with(cor(IC50_Coverage, SCC)) %>% round(3) %>% print   # PCC=-0.003
  IC50_Stat_Cell %>% plot_def(IC50_Coverage, SCC, main=main[1], xlab=xlab, ylab=ylab_c, add=add,
                              size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)
  
  
  ### 4. IC50 Mean
  print("# 4. IC50 Mean")
  ylab_c = bquote(SCC[C])
  ylab_d = bquote(SCC[D])
  xlab = bquote(Mean~ln(IC[50]))
  main = sprintf("%s/IC50_Mean & SCC Relation [%s]", dir, by)
  
  IC50_Stat_Drug %>% with(cor(IC50_Mean, SCC)) %>% round(3) %>% print   # PCC=-0.305
  IC50_Stat_Drug %>% plot_def(IC50_Mean, SCC, main=main[2], xlab=xlab, ylab=ylab_d, add=add,
                              size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)
  
  IC50_Stat_Cell %>% with(cor(IC50_Mean, SCC)) %>% round(3) %>% print   # PCC=0.220
  IC50_Stat_Cell %>% plot_def(IC50_Mean, SCC, main=main[1], xlab=xlab, ylab=ylab_c, add=add,
                              size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)
  
  
  ### 5. Num of Cell & Drug
  print("# 5. Num of Cell & Drug")
  subtitle = c("Cell", "Drug", "Cell TCGA", "Drug Pathway")
  main = sprintf("%s/IC50 Number & SCC [%s]", dir, subtitle)
  
  xlab = bquote(Number~of~ln(IC[50]))
  ylab_c = bquote(SCC[C])
  ylab_d = bquote(SCC[D])
  ylab_c_ = "SCC per Cancer Type"
  ylab_d_ = "SCC per Target Pathway"
  
  IC50_Stat_Drug %>% with(cor(Num, SCC)) %>% round(3) %>% print   # PCC=-0.115
  IC50_Stat_Drug %>% plot_def(Num, SCC, main=main[2], xlab=xlab, ylab=ylab_d,
                              size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)
  
  IC50_Stat_Cell %>% with(cor(Num, SCC)) %>% round(3) %>% print   # PCC=-0.198
  IC50_Stat_Cell %>% plot_def(Num, SCC, main=main[1], xlab=xlab, ylab=ylab_c,
                              size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)
  
  
  ### 6. Num of Cell & Drug [Group Level]
  print("# 6. Num of Cell & Drug [Group Level]")
  add = list(theme(axis.title.x=element_text(size=30)))
  
  IC50_Stat_Drug_Group %>% with(cor(Num, SCC)) %>% round(3) %>% print   # PCC=0.067
  IC50_Stat_Drug_Group %>% plot_def(Num, SCC, main=main[4], xlab=xlab, ylab=ylab_d_,
                                    size=2, alpha=0.5, axis_tl=25, axis_tx=18, add=add, save=T)
  
  IC50_Stat_Cell_Group %>% with(cor(Num, SCC)) %>% round(3) %>% print   # PCC=-0.11
  IC50_Stat_Cell_Group %>% plot_def(Num, SCC, main=main[3], xlab=xlab, ylab=ylab_c_,
                                    size=2, alpha=0.5, axis_tl=25, axis_tx=18, add=add, save=T)
  
  
  ### 7. IC50 Variation [GDSC-CCLE]
  print("# 7. IC50 Variation [GDSC-CCLE]")
  col = c("Cell_SANGER_ID", "Drug_CID", "LN_IC50_CCLE")
  by = c("Cell"="Cell_SANGER_ID", "Drug"="Drug_CID")
  Pred_GCNPath_CCLE = Pred_GCNPath_ %>% inner_join(IC50_CCLE[, col], by=by)   # 5567
  Stat_CCLE_List[[ttype]] = Pred_GCNPath_CCLE
  
  # min_ccle = min(Pred_GCNPath_CCLE$LN_IC50_CCLE)
  # max_ccle = max(Pred_GCNPath_CCLE$LN_IC50_CCLE)
  # sum(Pred_GCNPath_CCLE$LN_IC50_CCLE>=max_ccle)   # 2967 [53.29%]
  # sum(Pred_GCNPath_CCLE$LN_IC50_CCLE<=min_ccle)   # 38 [0.68%]
  #
  # Pred_GCNPath_CCLE$Cell %>% unique %>% length    # 374
  # Pred_GCNPath_CCLE$Drug %>% unique %>% length    # 18
  #
  # Pred_GCNPath_CCLE$LN_IC50 %>% is.na %>% sum         # 0
  # Pred_GCNPath_CCLE$LN_IC50_CCLE %>% is.na %>% sum    # 0
  # Pred_GCNPath_CCLE$LN_IC50_CCLE %>% hist             # Almost 2.07944154 [Max 8uM]
  
  
  xlab = bquote(GDSC~ln(IC[50]))
  ylab = bquote(CCLE~ln(IC[50]))
  main = sprintf("%s/LN_IC50 [GDSC-CCLE] & Pred_Error (color)", dir)
  
  add = list(scale_color_gradient(low="beige", high="firebrick3"))
  # midpoint = Pred_GCNPath_CCLE %>% with(abs(Prediction-LN_IC50)) %>% quantile(0.01)
  
  Pred_GCNPath_CCLE %>%
    mutate(Pred_Error=abs(Prediction-LN_IC50)) %>%
    plot_def(LN_IC50, LN_IC50_CCLE, color=Pred_Error,
             main=main, xlab=xlab, ylab=ylab, size=2.5, alpha=0.8,
             legend="|Error|", color_line="red", add=add,
             axis_tl=30, axis_tx=24, legend_tl=20, legend_tx=20, margin_lg=0.4,
             width=20, height=16.5, dpi=1000, xy_line=T, raster=T, save=T, save_svg=T)
  
  xlab = bquote("|"~GDSC~ln(IC[50]) - " " * CCLE~ln(IC[50])~"|")
  ylab = bquote("|"~GDSC~ln(IC[50]) - " " * Prediction~"|")
  main = sprintf("%s/LN_IC50 [GDSC-CCLE] & Pred_Error (no color)", dir)
  
  Pred_GCNPath_CCLE %>%
    mutate(Pred_Error=abs(Prediction-LN_IC50),
           LN_IC50_Var=abs(LN_IC50-LN_IC50_CCLE)) %>%
    plot_def(LN_IC50_Var, Pred_Error, main=main, xlab=xlab, ylab=ylab,
             size=2.5, alpha=0.5, axis_tl=22.5, axis_tx=22.5,
             width=16, height=16, dpi=1000, raster=T, save=T, save_svg=T)
  
  Pred_GCNPath_CCLE %>%
    mutate(Pred_Error=abs(Prediction-LN_IC50),
           LN_IC50_Var=abs(LN_IC50-LN_IC50_CCLE)) %>%
    with(cor(Pred_Error, LN_IC50_Var)) %>% round(3) %>% print
}

by = c("Cell", "Drug")
xlab = bquote(Mean~ln(IC[50]))
ylab = bquote(SD~ln(IC[50]))
dir = mkdir("Analysis of Error")
main = sprintf("%s/IC50_Mean & IC50_SD Relation [GDSC, %s]", dir, by)

col = c("Drug", "Drug_Pathway", "IC50_Mean", "IC50_SD", "Num")
IC50_Stat_Drug_ = IC50_Stat_Drug[, col]
IC50_Stat_Drug_ %>% with(cor(IC50_Mean, IC50_SD)) %>% round(3)   # -0.45
IC50_Stat_Drug_ %>% plot_def(IC50_Mean, IC50_SD, main=main[2], xlab=xlab, ylab=ylab, add=add,
                             size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)

col = c("Cell", "Cell_TCGA", "IC50_Mean", "IC50_SD", "Num")
IC50_Stat_Cell_ = IC50_Stat_Cell[, col]
IC50_Stat_Cell_ %>% with(cor(IC50_Mean, IC50_SD)) %>% round(3)   # -0.407
IC50_Stat_Cell_ %>% plot_def(IC50_Mean, IC50_SD, main=main[1], xlab=xlab, ylab=ylab, add=add,
                             size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)

IC50_Stat_Drug_CCLE = IC50_CCLE %>% group_by(Drug) %>% 
  summarise(Num=n(), IC50_Mean=mean(LN_IC50_CCLE), 
            IC50_SD=sd(LN_IC50_CCLE)) %>% as.data.frame   # 24

IC50_Stat_Cell_CCLE = IC50_CCLE %>% group_by(Cell) %>% 
  summarise(Num=n(), IC50_Mean=mean(LN_IC50_CCLE), 
            IC50_SD=sd(LN_IC50_CCLE)) %>% as.data.frame   # 504

main = sprintf("%s/IC50_Mean & IC50_SD Relation [CCLE, %s]", dir, by)
IC50_Stat_Drug_CCLE %>% with(cor(IC50_Mean, IC50_SD)) %>% round(3)   # -0.22
IC50_Stat_Drug_CCLE %>% 
  plot_def(IC50_Mean, IC50_SD, main=main[2], xlab=xlab, ylab=ylab, add=add,
           size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)

IC50_Stat_Cell_CCLE %>% with(cor(IC50_Mean, IC50_SD)) %>% round(3)   # -0.241
IC50_Stat_Cell_CCLE %>% 
  plot_def(IC50_Mean, IC50_SD, main=main[1], xlab=xlab, ylab=ylab, add=add,
           size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)



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
  
  process_perf = function(Perf_List, score="RMSE") {
    Perf_List_ = Perf_List %>% 
      rename(Train_Fold=Fold, Num_Test=N_Test) %>% 
      mutate(Dataset=recode(Dataset, "GDSC"="GDSC1+2")) %>% 
      subset(select=c(Model, Dataset, Test_Type, Train_Fold, Num_Test, object(score)))   # 2310
    
    test_type = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")
    Perf_List_N = Perf_List_ %>% subset(Test_Type==test_type[1]) %>% subset(select=-Test_Type)
    Perf_List_C = Perf_List_ %>% subset(Test_Type==test_type[2]) %>% subset(select=-Test_Type)
    Perf_List_D = Perf_List_ %>% subset(Test_Type==test_type[3]) %>% subset(select=-Test_Type)
    Perf_List_S = Perf_List_ %>% subset(Test_Type==test_type[4]) %>% subset(select=-Test_Type)
    
    Perf_List_ = list(Perf_List_N, Perf_List_C, Perf_List_D, Perf_List_S)
    return(Perf_List_)
  }
  
  process_utest = function(Perf_UTest, score="RMSE") {
    score_ = sprintf("_%s", score)
    col_rm = c("RMSE", "PCC", "SCC") %>% setdiff(score)
    
    Perf_UTest_ = Perf_UTest %>% 
      subset(Dataset=="GDSC") %>% subset(select=-Dataset) %>% 
      select(!contains(col_rm)) %>% rename_with(~gsub(score_, "", .x)) %>% 
      rename(Minus_Log10_FDR=MLog10_FDR) %>% 
      relocate(Minus_Log10_FDR, .after=FDR) %>% as.data.frame
    
    test_type = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")
    Perf_UTest_N = Perf_UTest_ %>% subset(Test_Type==test_type[1]) %>% subset(select=-Test_Type)
    Perf_UTest_C = Perf_UTest_ %>% subset(Test_Type==test_type[2]) %>% subset(select=-Test_Type)
    Perf_UTest_D = Perf_UTest_ %>% subset(Test_Type==test_type[3]) %>% subset(select=-Test_Type)
    Perf_UTest_S = Perf_UTest_ %>% subset(Test_Type==test_type[4]) %>% subset(select=-Test_Type)
    
    Perf_UTest_List = list(Perf_UTest_N, Perf_UTest_C, Perf_UTest_D, Perf_UTest_S)
    return(Perf_UTest_List)
  }
  
  col = c("MODEL_NAME", "SANGER_MODEL_ID")
  col_ = c("Cell_Name", "Cell_SANGER_ID")
  Anno_Cells_ = Anno_Cells[, col] %>% setNames(col_) %>% 
    subset(!is.na(Cell_SANGER_ID)) %>% distinct(Cell_SANGER_ID, .keep_all=T) %>% as.data.frame
  
  col = c("Name", "Drug_CID")
  col_ = c("Drug_Name", "Drug_CID")
  Anno_Drugs_ = Anno_Drugs[, col] %>% setNames(col_) %>% 
    subset(!is.na(Drug_CID)) %>% distinct(Drug_CID, .keep_all=T) %>% as.data.frame
  
  test_type_ = c("Unblinded", "Cell_Blind", "Drug_Blind", "Strict_Blind")
  
  ### [Source Data] Supplementary Fig. 4-7
  if (save_pred) {
    Pred_U %>% save_for_nc(num=4, suppl=T)
    Pred_C %>% save_for_nc(num=5, suppl=T)
    Pred_D %>% save_for_nc(num=6, suppl=T)
    Pred_S %>% save_for_nc(num=7, suppl=T)
  }
  
  ### [Source Data] Supplementary Fig. 8
  Perf_RMSE_ = Perf_List %>% process_perf(score="RMSE")
  Perf_RMSE_ %>% save_for_nc(num=8, suppl=T)
  
  ### [Source Data] Supplementary Fig. 9
  Perf_PCC_ = Perf_List %>% process_perf(score="PCC")
  Perf_PCC_ %>% save_for_nc(num=9, suppl=T)
  
  ### [Source Data] Supplementary Fig. 10
  Perf_SCC_ = Perf_List %>% process_perf(score="SCC")
  Perf_SCC_ %>% save_for_nc(num=10, suppl=T)
  
  ### [Source Data] Supplementary Fig. 11
  Perf_UTest_RMSE_ = Perf_UTest %>% process_utest(score="RMSE")
  Perf_UTest_RMSE_ %>% save_for_nc(num=11, suppl=T)
  
  ### [Source Data] Supplementary Fig. 12
  Perf_UTest_PCC_ = Perf_UTest %>% process_utest(score="PCC")
  Perf_UTest_PCC_ %>% save_for_nc(num=12, suppl=T)
  
  ### [Source Data] Supplementary Fig. 13
  Perf_UTest_SCC_ = Perf_UTest %>% process_utest(score="SCC")
  Perf_UTest_SCC_ %>% save_for_nc(num=13, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 14
  Perf_Cell_TCGA_ = Perf_Cell_TCGA %>% 
    subset(Dataset=="GDSC" & Test_Type %in% test_type[c(1, 2, 4)]) %>% 
    mutate(Test_Type=ifelse(Test_Type!="Normal", Test_Type, "Unblinded")) %>% 
    rename(TCGA_Code=Cell_TCGA) %>% subset(select=-Dataset)
  
  Perf_Cell_TCGA_ %>% save_for_nc(num=14, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 15
  Perf_Drug_Pathway_ = Perf_Drug_Pathway %>% 
    subset(Dataset=="GDSC" & Test_Type %in% test_type[c(1, 3, 4)]) %>% 
    mutate(Test_Type=ifelse(Test_Type!="Normal", Test_Type, "Unblinded")) %>% 
    rename(Target_Pathway=Drug_Pathway) %>% subset(select=-Dataset)
  
  Perf_Drug_Pathway_ %>% save_for_nc(num=15, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 16
  Perf_Cell_ = Perf_Cell %>% 
    subset(Model=="GCNPath" & Dataset=="GDSC" & Test_Type %in% test_type[c(1, 2, 4)]) %>% 
    mutate(Test_Type=ifelse(Test_Type!="Normal", Test_Type, "Unblinded") %>% factor(levels=test_type_)) %>% 
    arrange(Test_Type, Cell) %>% rename(Num_Test=N_Test) %>% subset(select=-c(Model, Dataset))
  
  by = c("Cell"="Cell_SANGER_ID")
  Perf_Cell_ = Perf_Cell_ %>% 
    left_join(Anno_Cells_, by=by) %>% 
    rename(Cell_SANGER_ID=Cell, TCGA_Code=Cell_TCGA) %>% 
    relocate(Cell_Name, .before=everything())
  
  Temp = list(Perf_Cell_ %>% subset(Test_Type==test_type_[1]) %>% subset(select=-Test_Type), 
              Perf_Cell_ %>% subset(Test_Type==test_type_[2]) %>% subset(select=-Test_Type), 
              Perf_Cell_ %>% subset(Test_Type==test_type_[4]) %>% subset(select=-Test_Type))
  
  Temp %>% save_for_nc(num=16, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 17
  Perf_Drug_ = Perf_Drug %>% 
    subset(Model=="GCNPath" & Dataset=="GDSC" & Test_Type %in% test_type[c(1, 3, 4)]) %>% 
    mutate(Test_Type=ifelse(Test_Type!="Normal", Test_Type, "Unblinded") %>% factor(levels=test_type_)) %>% 
    arrange(Test_Type, Drug) %>% rename(Num_Test=N_Test) %>% subset(select=-c(Model, Dataset))
  
  by = c("Drug"="Drug_CID")
  Perf_Drug_ = Perf_Drug_ %>% 
    left_join(Anno_Drugs_, by=by) %>% 
    rename(Drug_CID=Drug, Target_Pathway=Drug_Pathway) %>% 
    relocate(Drug_Name, .before=everything())
  
  Temp = list(Perf_Drug_ %>% subset(Test_Type==test_type_[1]) %>% subset(select=-Test_Type), 
              Perf_Drug_ %>% subset(Test_Type==test_type_[3]) %>% subset(select=-Test_Type), 
              Perf_Drug_ %>% subset(Test_Type==test_type_[4]) %>% subset(select=-Test_Type))
  
  Temp %>% save_for_nc(num=17, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 18
  Perf_Cell_TCGA_ = Perf_Cell_TCGA %>% 
    subset(Dataset=="GDSC1" & Test_Type %in% test_type[c(1, 2, 4)]) %>% 
    mutate(Test_Type=ifelse(Test_Type!="Normal", Test_Type, "Unblinded")) %>% 
    rename(TCGA_Code=Cell_TCGA) %>% subset(select=-Dataset)
  
  Perf_Drug_Pathway_ = Perf_Drug_Pathway %>% 
    subset(Dataset=="GDSC1" & Test_Type %in% test_type[c(1, 3, 4)]) %>% 
    mutate(Test_Type=ifelse(Test_Type!="Normal", Test_Type, "Unblinded")) %>% 
    rename(Target_Pathway=Drug_Pathway) %>% subset(select=-Dataset)
  
  Perf_GCNPath_GDSC1C_ = Perf_GCNPath_GDSC1C %>% rename(TCGA_Code=Cell_TCGA)
  Perf_GCNPath_GDSC1D_ = Perf_GCNPath_GDSC1D %>% rename(Target_Pathway=Drug_Pathway)
  
  Temp = list(Perf_Cell_TCGA_, Perf_GCNPath_GDSC1C_, 
              Perf_Drug_Pathway_, Perf_GCNPath_GDSC1D_)
  
  num_fig = c("a_top", "a_bottom", "b_top", "b_bottom")
  Temp %>% save_for_nc(num=18, num_fig=num_fig, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 19
  Perf_Cell_TCGA_ = Perf_Cell_TCGA %>% 
    subset(Dataset=="GDSC2" & Test_Type %in% test_type[c(1, 2, 4)]) %>% 
    mutate(Test_Type=ifelse(Test_Type!="Normal", Test_Type, "Unblinded")) %>% 
    rename(TCGA_Code=Cell_TCGA) %>% subset(select=-Dataset)
  
  Perf_Drug_Pathway_ = Perf_Drug_Pathway %>% 
    subset(Dataset=="GDSC2" & Test_Type %in% test_type[c(1, 3, 4)]) %>% 
    mutate(Test_Type=ifelse(Test_Type!="Normal", Test_Type, "Unblinded")) %>% 
    rename(Target_Pathway=Drug_Pathway) %>% subset(select=-Dataset)
  
  Perf_GCNPath_GDSC2C_ = Perf_GCNPath_GDSC2C %>% rename(TCGA_Code=Cell_TCGA)
  Perf_GCNPath_GDSC2D_ = Perf_GCNPath_GDSC2D %>% rename(Target_Pathway=Drug_Pathway)
  
  Temp = list(Perf_Cell_TCGA_, Perf_GCNPath_GDSC2C_, 
              Perf_Drug_Pathway_, Perf_GCNPath_GDSC2D_)
  
  num_fig = c("a_top", "a_bottom", "b_top", "b_bottom")
  Temp %>% save_for_nc(num=19, num_fig=num_fig, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 23
  Perf_Cell_Corr_ = Perf_Cell_Corr %>% 
    lapply(function(df) df[, c(1, 2, 5)] %>% 
             rename(Model_X=Model1, Model_Y=Model2, PCC_of_SCC=SCC_Corr))
  Perf_Cell_Corr_ %>% save_for_nc(num=23, suppl=T)
  
  ### [Source Data] Supplementary Fig. 24
  Perf_Drug_Corr_ = Perf_Drug_Corr %>% 
    lapply(function(df) df[, c(1, 2, 5)] %>% 
             rename(Model_X=Model1, Model_Y=Model2, PCC_of_SCC=SCC_Corr))
  Perf_Drug_Corr_ %>% save_for_nc(num=24, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 25
  Stat_Cell_List_ = Stat_Cell_List[c(1, 2, 4)] %>%
    lapply(function(df) df %>%
             dplyr::rename(TCGA_Code=Cell_TCGA, Num_Test=Num,
                           LN_IC50_Mean=IC50_Mean, LN_IC50_SD=IC50_SD,
                           LN_IC50_Bimodality=IC50_Bimodality,
                           LN_IC50_Coverage=IC50_Coverage))
  
  Stat_Drug_List_ = Stat_Drug_List[c(1, 3, 4)] %>%
    lapply(function(df) df %>%
             dplyr::rename(Target_Pathway=Drug_Pathway, Num_Test=Num,
                           LN_IC50_Mean=IC50_Mean, LN_IC50_SD=IC50_SD,
                           LN_IC50_Bimodality=IC50_Bimodality,
                           LN_IC50_Coverage=IC50_Coverage))
  
  Stat_Cell_List_Group_ = Stat_Cell_List_Group[c(1, 2, 4)] %>%
    lapply(function(df) df %>% dplyr::rename(TCGA_Code=Cell_TCGA, Num_Test=Num))
  
  Stat_Drug_List_Group_ = Stat_Drug_List_Group[c(1, 3, 4)] %>%
    lapply(function(df) df %>% dplyr::rename(Target_Pathway=Drug_Pathway, Num_Test=Num))
  
  IC50_Stat_ = list(Stat_Drug_List_[[1]], Stat_Drug_List_Group_[[1]], 
                    Stat_Drug_List_[[2]], Stat_Drug_List_Group_[[2]], 
                    Stat_Drug_List_[[3]], Stat_Drug_List_Group_[[3]], 
                    Stat_Cell_List_[[1]], Stat_Cell_List_Group_[[1]], 
                    Stat_Cell_List_[[2]], Stat_Cell_List_Group_[[2]],
                    Stat_Cell_List_[[3]], Stat_Cell_List_Group_[[3]])
  
  direction = c("left", "right")
  num_fig = expand.grid(letters[1:6], direction) %>% arrange(Var1)
  num_fig = num_fig %>% with(sprintf("%s_%s", Var1, Var2))
  IC50_Stat_ %>% save_for_nc(num=25, num_fig=num_fig, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 26
  Stat_CCLE_List_ = Stat_CCLE_List %>% 
    lapply(function(df) df %>% 
             subset(select=-c(Model, Dataset, Test_Type, Cell_BROAD, Cell_COSMIC, Cell_TCGA, Drug_Pathway)) %>% 
             rename(Train_Fold=Fold, LN_IC50_GDSC=LN_IC50) %>% 
             mutate(Diff_GDSC_Pred=abs(Prediction-LN_IC50_GDSC), 
                    Diff_GDSC_CCLE=abs(LN_IC50_CCLE-LN_IC50_GDSC)) %>% 
             relocate(Train_Fold, .after=Drug) %>% 
             relocate(Prediction, .after=LN_IC50_CCLE) %>% as.data.frame)
  
  Stat_CCLE_List_ %>% save_for_nc(num=26, suppl=T)
  
  
  ### [Source Data] Fig. 3
  subset_db = function(Perf) Perf %>% subset(Dataset=="GDSC1+2") %>% subset(select=-Dataset)
  Perf_RMSE_GDSC_ = Perf_RMSE_ %>% lapply(subset_db)
  Perf_RMSE_GDSC_ %>% sapply(nrow)   # 140, 140, 140, 350
  Perf_RMSE_GDSC_ %>% save_for_nc(num=3, suppl=F)
}
