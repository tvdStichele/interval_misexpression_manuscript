#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/sv_carrier_count_gene_body/carrier_info_chr%I.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/sv_carrier_count_gene_body/carrier_info_chr%I.out
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
paired_smpls=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/1_rna_seq_qc/wgs_rna_match/paired_wgs_rna_postqc_prioritise_wgs.tsv
express=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/2_misexp_qc/misexp_gene_cov_corr/tpm_zscore_4568smpls_8610genes_flat_misexp_corr_qc.csv
vep_msc=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/4_vrnt_enrich/sv_vep/msc/SV_vep_hg38_msc_parsed.tsv
gencode=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/reference/gencode/gencode.v31.annotation.sorted.gtf.gz
sv_info=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/data/sv_vcf/info_table/final_sites_critical_info_allele.txt
root=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/4_vrnt_enrich/sv_count_carriers/gene_body/200kb_window

python ./sv_carrier_express.py --express $express \
                               --chrom chr$CHR \
                               --sv_info $sv_info \
                               --vcf $vcf \
                               --gencode $gencode \
                               --paired_smpls $paired_smpls \
                               --vep_msc $vep_msc \
                               --window 200000 --root $root