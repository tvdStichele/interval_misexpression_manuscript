#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/vrnts_bed/vrnts2bed_%I.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/vrnts_bed/vrnts2bed_%I.out
#BSUB -q "basement"
#BSUB -R "select[mem>50000] rusage[mem=50000] span[hosts=1]"
#BSUB -M50000
#BSUB -n 2
#BSUB -J vrnts2bed[1-22]

# conda activate tv5_base
i=`expr ${LSB_JOBINDEX} - 1`
chroms=({1..22})
CHR=${chroms[$i]}

vcf=/lustre/scratch126/humgen/projects/interval_wgs/final_release_freeze/gt_phased/interval_wgs.chr${CHR}.gt_phased.vcf.gz
root=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/4_vrnt_enrich/snp_indel_count_carriers

python vrnts2bed.py --chrom chr$CHR --vcf $vcf --max_indel_length 50 --root $root