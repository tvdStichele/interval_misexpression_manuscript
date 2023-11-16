#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/bam_extract/vrnt%I_misexp_cntrl_bams.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/bam_extract/vrnt%I_misexp_cntrl_bams.out
#BSUB -q "basement"
#BSUB -R "select[mem>50000] rusage[mem=50000] span[hosts=1]"
#BSUB -M50000
#BSUB -n 1
#BSUB -J extract_bam[1-105]

# conda activate tv5_base
wkdir=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3

star_dir=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq_n5188_v97/results/star_pass2_2ndpass/
misexp_vrnts=${wkdir}/5_misexp_vrnts/test_cntrl_sets/vrnt_id_misexp_tpm_zscore_median.txt
misexp_feat=${wkdir}/6_misexp_dissect/vrnt_features/misexp_vrnt_features.tsv
rna_id_pass=${wkdir}/1_rna_seq_qc/aberrant_smpl_qc/smpls_pass_qc_4568.csv
out=${wkdir}/6_misexp_dissect/bams
mkdir -p $out

## get variant from list
awk 'NR=='${LSB_JOBINDEX} $misexp_vrnts | while read -a line
do
    vrnt_id=${line[0]}
    python bam_file_extract.py --vrnt_id $vrnt_id --star_dir $star_dir --rna_pass $rna_id_pass --misexp_feat $misexp_feat --out $out
done
