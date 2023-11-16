#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/cadd_sv/cadd_sv_intrvl_all_svs_%I.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/cadd_sv/cadd_sv_intrvl_all_svs_%I.out
#BSUB -q "basement"
#BSUB -R "select[mem>10000] rusage[mem=10000] span[hosts=1]"
#BSUB -M10000
#BSUB -n 4
#BSUB -J cadd_sv_intrvl_all[1-25]

# conda activate run.caddsv
wkdir=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3
configs=${wkdir}/5_misexp_vrnts/scores/cadd_sv/config_paths_no_invs.txt
awk 'NR=='${LSB_JOBINDEX} $configs | while read -a line
do 
    config_path=${line[0]}
    cd /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/software/CADD-SV
    snakemake  --use-conda --configfile $config_path -j 4
done 