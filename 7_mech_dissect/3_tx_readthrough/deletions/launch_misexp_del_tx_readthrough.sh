#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/del_tx_read_intergenic/vrnt%I_tx_read_intergenic.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/del_tx_read_intergenic/vrnt%I_tx_read_intergenic.out
#BSUB -q "normal"
#BSUB -R "select[mem>40000] rusage[mem=40000] span[hosts=1]"
#BSUB -M40000
#BSUB -n 1
#BSUB -J del_tx_read[1-12]

# conda activate tv5_base
wkdir=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3

tx_read_candidates=${wkdir}/6_misexp_dissect/tx_readthrough/deletions/misexp_del_tx_vrnts.txt
tx_read_info_path=${wkdir}/6_misexp_dissect/tx_readthrough/deletions/misexp_del_tx_readthrough_candidates.tsv
bams_dir=${wkdir}/6_misexp_dissect/bams
out=${wkdir}/6_misexp_dissect/tx_readthrough/deletions/readthrough_region

## get variant from list
awk 'NR=='${LSB_JOBINDEX} $tx_read_candidates | while read -a line
do
    vrnt_id=${line[0]}
    python misexp_del_tx_read_intergenic.py --vrnt_id $vrnt_id --tx_read_info $tx_read_info_path --bam_dir $bams_dir --out $out
done
