#!/bin/bash
#BSUB -o /lustre/scratch119/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression/4_vrnt_enrich/1_snp_indel/flat_vcf/logs/vrnts2flatvcf_%I.out
#BSUB -e /lustre/scratch119/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression/4_vrnt_enrich/1_snp_indel/flat_vcf/logs/vrnts2flatvcf_%I.out
#BSUB -q "normal"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -M100000
#BSUB -n 4
#BSUB -J vrnts2bed[21]

# conda activate /lustre/scratch119/realdata/mdt3/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/envs/misexp_py
i=`expr ${LSB_JOBINDEX} - 1`
chroms=({1..22})
CHR=${chroms[$i]}

vcf=/lustre/scratch119/humgen/projects/interval_wgs/release/chr${CHR}/chr${CHR}.intervalwgs_v2_GT_only.eagle_phased.vcf.gz
root=/lustre/scratch119/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression/4_vrnt_enrich/1_snp_indel/

python vrnts2flatvcf.py --chrom chr$CHR \
                        --vcf $vcf \
                        --root $root