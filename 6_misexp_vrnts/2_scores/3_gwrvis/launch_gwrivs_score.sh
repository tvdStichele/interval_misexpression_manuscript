#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/gwrvis/gwrvis_score_chr%I.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/gwrvis/gwrvis_score_chr%I.out
#BSUB -q "normal"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -M100000
#BSUB -n 1
#BSUB -J gwrvis_score[1-22]

# conda activate tv5_base

i=`expr ${LSB_JOBINDEX} - 1`
chroms=({1..22})
CHR=${chroms[$i]}

wkdir=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3
gwrvis=/lustre/scratch125/humgen/resources/JARVIS-gwRVIS-scores/hg38/gwRVIS/gwrvis_single_nt.chr${CHR}.hg38.bed.gz
sv_info=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/data/sv_vcf/info_table/final_sites_critical_info_allele.txt
sv_bed_path=${wkdir}/5_misexp_vrnts/test_cntrl_sets/vrnt_id_in_windows_chr_num_misexp_genes.bed
out=${wkdir}/5_misexp_vrnts/scores/gwrvis

python sv_gwrvis_score.py --chrom chr$CHR --gwrvis $gwrvis --sv_info $sv_info --sv_bed $sv_bed_path --out $out