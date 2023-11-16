#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/logs/sv_carrier_windows/sv_carrier_count_chr%I.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/logs/sv_carrier_windows/sv_carrier_count_chr%I.out
#BSUB -q "normal"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -M100000
#BSUB -n 1
#BSUB -J sv_carrier_count[1-22]

# conda activate tv5_base
i=`expr ${LSB_JOBINDEX} - 1`
chroms=({1..22})
CHR=${chroms[$i]}

vcf=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/data/sv_vcf/filtered_merged_gs_svp_10728.vcf.gz
paired_smpls=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/1_rna_seq_qc/wgs_rna_match/paired_wgs_rna_postqc_prioritise_wgs.tsv
express=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/2_misexp_qc/misexp_gene_cov_corr/tpm_zscore_4568smpls_8610genes_flat_misexp_corr_qc.csv
gencode=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/reference/gencode/gencode.v31.annotation.sorted.gtf.gz
sv_info=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/data/sv_vcf/info_table/final_sites_critical_info_allele.txt
vep_msc=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/4_vrnt_enrich/sv_vep/msc/SV_vep_hg38_msc_parsed.tsv
root=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/4_vrnt_enrich/sv_count_carriers/windows/200kb_window

python ./count_sv_carrier_window_zscore_msc.py --express $express \
                                   --chrom chr$CHR \
                                   --sv_info $sv_info \
                                   --vcf $vcf \
                                   --gencode $gencode \
                                   --paired_smpls $paired_smpls \
                                   --z_cutoff 2 5 10 15 20 25 30 35 40 45 50 --window_start 0 \
                                   --window_step_size 200000 \
                                   --window_max 1000000 \
                                   --vep_msc $vep_msc \
                                   --af_range 0 0.01 --root $root