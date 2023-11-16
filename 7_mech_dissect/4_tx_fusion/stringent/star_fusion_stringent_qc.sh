#!/bin/bash
sample=$1
genome_lib_dir=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression/8_vignette/star_fusion/ctat_genome_lib/output
left_fq="/lustre/scratch126/humgen/projects/interval_rna/nf_ci_cram_fastqs_5591/results/crams_to_fastq/fastq/5591/${sample}/5591.${sample}_1.fastq.gz"
right_fq="/lustre/scratch126/humgen/projects/interval_rna/nf_ci_cram_fastqs_5591/results/crams_to_fastq/fastq/5591/${sample}/5591.${sample}_2.fastq.gz"
output="/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/6_misexp_dissect/star_fusion/results_stringent_qc/${sample}"
mkdir -p $output

/usr/local/src/STAR-Fusion/STAR-Fusion --genome_lib_dir $genome_lib_dir \
                            --left_fq $left_fq \
                            --right_fq $right_fq \
                            --output_dir $output \
                            --STAR_max_mate_dist 50000 \
                            --FusionInspector validate \
                            --denovo_reconstruct \
                            --examine_coding_effect \
                            --no_annotation_filter # disable annotation-based filtering 
                            
             