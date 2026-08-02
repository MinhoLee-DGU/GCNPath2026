#!/usr/bin/bash

gpu=0
node=lys
test_type=Normal

ic50=data/ic50_data/IC50_SCLC.txt
cell=processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle
drug=processed/drug_data/SCLC_Drug_Graph.pickle
dir_model=results/IC50_GDSC/$test_type/RGCN

col_cell=Cell
col_drug=Drug
col_ic50=0

use_slurm=0
seed_list=$(seq 2021 2030)
for seed in ${seed_list[@]}
do
    jname=SCLC_Seed${seed}
    param=${dir_model}/param_retrain_seed${seed}.pt
    hparam=${dir_model}/hyper_param_retrain_seed${seed}.pickle
    out_file=${dir_model}/pred_sclc_seed${seed}.csv
    out_cam=${dir_model}/pred_sclc_seed${seed}_gcam.csv
    bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
        "$out_cam" $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
done

