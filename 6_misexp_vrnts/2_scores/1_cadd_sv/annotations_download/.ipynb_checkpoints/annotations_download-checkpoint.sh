#!/bin/bash
#BSUB -o /lustre/scratch119/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/phenotypes/rna_seq/tests/derepression/cadd_sv/annotations/logs/annotations.out
#BSUB -e /lustre/scratch119/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/phenotypes/rna_seq/tests/derepression/cadd_sv/annotations/logs/annotations.err
#BSUB -q "long"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -M100000
#BSUB -n 4
#BSUB -J cadd_sv_annotations

cd /lustre/scratch119/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/software/CADD-SV
# download annotation folders 
wget https://kircherlab.bihealth.org/download/CADD-SV/v1.0/dependencies.tar.gz