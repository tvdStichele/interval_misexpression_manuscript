#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/sv_count_carriers_10kb/sv_carrier_count_gene_msc_reg_chr%I.out
#BSUB -e  /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/sv_count_carriers_10kb/sv_carrier_count_gene_msc_reg_chr%I.out
#BSUB -q "normal"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -M100000
#BSUB -n 1
#BSUB -J sv_carrier_count[1-22]

# conda activate tv5_base
i=`expr ${LSB_JOBINDEX} - 1`
chroms=({1..22})
CHR=${chroms[$i]}

vep_msc=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/4_vrnt_enrich/sv_vep/msc/SV_vep_hg38_msc_parsed.tsv
vep_all=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/4_vrnt_enrich/sv_vep/all/SV_vep_hg38_all_parsed.tsv
sv_info=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/data/sv_vcf/info_table/final_sites_critical_info_allele.txt
root=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/4_vrnt_enrich/sv_count_carriers/gene_body/10kb_window

python ./count_sv_carrier_vrnt_gene_msc_reg.py --chrom chr$CHR \
                                               --sv_info $sv_info \
                                               --vep_all $vep_all \
                                               --vep_msc $vep_msc \
                                               --z_cutoff 2 3 4 5 10 15 20 25 30 35 40 45 50\
                                               --af_range 0 0.01 0.05 0.1 0.5 --root $root