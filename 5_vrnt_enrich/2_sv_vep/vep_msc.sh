#!/usr/bin/env bash
# inputs
VCF_IN=$1
OUT_DIR=$2

VCF_FILE=$(basename "$VCF_IN")
VCF_DIR=$(dirname "$VCF_IN")

set -o errexit
set -o pipefail
set -o nounset

oufnprfx="SV_vep_hg38_msc"
mkdir -p $OUT_DIR
cd $OUT_DIR
echo $VCF_DIR
echo $VCF_FILE
echo $VCF_IN

/software/singularity-v3.6.4/bin/singularity exec \
--bind $VCF_DIR:/opt/vcf \
--bind /lustre/scratch125/humgen/resources/ensembl/vep/GRCh38/vep_data:/opt/vep/.vep \
--bind /lustre/scratch125/humgen/resources/ensembl/vep/GRCh38/Plugins:/opt/vep/.vep/Plugins \
--bind /lustre/scratch125/humgen/resources/gnomAD/release-2.1.1/exomes \
--bind /lustre/scratch125/humgen/resources/cadd_scores/20201027-GRCh38_v1.6 \
--bind /lustre/scratch125/humgen/resources/SpliceAI_data_files \
--bind $OUT_DIR:/opt/vcf/out \
/lustre/scratch125/humgen/resources/ensembl/vep/singularity_containers/vep_97.3.sif \
vep \
-i /opt/vcf/$VCF_FILE \
-o /opt/vcf/out/${oufnprfx}_output.txt \
--offline \
--cache \
--dir_cache /opt/vep/.vep/ \
--verbose \
--assembly GRCh38 \
--fork 4 \
--most_severe \
--fasta /opt/vep/.vep/homo_sapiens/97_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --regulatory 
