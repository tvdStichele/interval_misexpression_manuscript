#!/bin/bash

set -euo pipefail

#User inputs
SCRIPT=$1
SAMPLE=$2
#All arguments from third onwards are put into the variable ARG
#ARGS=${@:3}

log_dir="/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/star_fusion_stringent_qc/${SAMPLE}"
mkdir -p $log_dir
#BSUB ARGUMENTS

CPUS=4
MEM=50000 # ~40G RAM required
QUE="normal"
GROUP="team282"
IMAGE="/software/CASM/singularity/star-fusion/star-fusion_1.10.1.sif"
LOG="${log_dir}/${SAMPLE}.out"

###################### DONT CHANGE OPTIONS BELOW THIS LINE ###########################

bsub \
  -G $GROUP \
  -n $CPUS \
  -R"span[hosts=1] select[mem>${MEM}] rusage[mem=${MEM}]" \
  -M $MEM \
  -o $LOG \
  -e $LOG \
  -q $QUE \
  -J $SAMPLE \
  /software/singularity-v3.9.0/bin/singularity exec -B /lustre,/nfs $IMAGE $SCRIPT $SAMPLE #$ARG