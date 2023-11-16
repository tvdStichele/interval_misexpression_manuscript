#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/annotate_vrnts_logs/chr%I_launch.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/annotate_vrnts_logs/chr%I_launch.out
#BSUB -q "normal"
#BSUB -R "select[mem>120000] rusage[mem=120000] span[hosts=1]"
#BSUB -M120000
#BSUB -n 1
#BSUB -J launch_chrom[3-7]

# conda activate tv5_base
i=`expr ${LSB_JOBINDEX} - 1`
chroms=({1..22})
CHR=${chroms[$i]}

script=annotate_vrnts_by_gene.sh
wkdir=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3

gene_id_list=${wkdir}/4_vrnt_enrich/snp_indel_count_carriers/intersect_genes/chr${CHR}_gene_ids_intersect.txt
num_jobs=`cat $gene_id_list | wc -l`

log_dir=${wkdir}/logs/annotate_by_gene_logs/chr${CHR}
mkdir -p $log_dir

bsub -o ${log_dir}/gene_%I.out -e ${log_dir}/gene_%I.out -q "normal" -R "select[mem>50000] rusage[mem=50000] span[hosts=1]" -M50000 -n 1 -J "annotate_vrnts[1-${num_jobs}]" "sh ${script} ${CHR} ${gene_id_list}"