#!/bin/bash
collapse_script=$1
out_dir=$2

### GENCODE 
gencode_http=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
gencode_dir=${out_dir}/gencode
mkdir -p $gencode_dir
gtf_file=$(basename $gencode_http)
gtf_name=$(echo $gtf_file | sed 's/\.gtf.gz$//')

# download, sort, bgzip and index gencode annotation file 
cd $gencode_dir
wget -nc $gencode_http
gtf_in=${gencode_dir}/${gtf_file}
gtf_sorted=${gencode_dir}/${gtf_name}.sorted.gtf
bedtools sort -i $gtf_in > $gtf_sorted
bgzip -f $gtf_sorted
tabix -fp gff ${gtf_sorted}.gz

# collapse gencode file using GTEx method, sort, bgzip and index output
gtf_collapse=${gencode_dir}/${gtf_name}.collapsed.gtf
python3 $collapse_script $gtf_in $gtf_collapse --collapse_only --stranded
gtf_collapse_sorted=${gencode_dir}/${gtf_name}.collapsed.sorted.gtf
bedtools sort -i $gtf_collapse > $gtf_collapse_sorted
rm $gtf_collapse
bgzip -f $gtf_collapse_sorted
tabix -fp gff ${gtf_collapse_sorted}.gz

### ENSEMBL
ensembl_http=http://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
ensembl_dir=${out_dir}/ensembl
mkdir -p $ensembl_dir
cd $ensembl_dir
ensembl_file=$(basename $ensembl_http)
ensembl_name=$(echo $ensembl_file | sed 's/\.gtf.gz$//')
# download, sort, bgzib and index ensembl gtf 
wget -nc $ensembl_http 
ensembl_in=${ensembl_dir}/${ensembl_file}
ensembl_sorted=${ensembl_dir}/${ensembl_name}.sorted.gtf
bedtools sort -i $ensembl_in > $ensembl_sorted
# bgzip 
bgzip -f $ensembl_sorted
# tabix file 
tabix -fp gff ${ensembl_sorted}.gz

### GTEX 
gtex_dir=${out_dir}/gtex
mkdir -p $gtex_dir

# download GTEx eQTL expression matrices data
gtex_eqtl_http=https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_expression_matrices.tar 
mkdir -p ${gtex_dir}/eqtl
cd ${gtex_dir}/eqtl
gtex_eqtl_filename=$(basename $gtex_eqtl_http)
wget -nc $gtex_eqtl_http
# untar and remove tar file
tar -xvf $gtex_eqtl_filename

# download GTEx European eQTL data 
gtex_eur_eqtl_http=https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_EUR.tar
gtex_eur_eqtl_filename=$(basename $gtex_eur_eqtl_http)
gtex_eur_eqtl_dir=$(echo $gtex_eur_eqtl_filename | sed 's/\.tar$//')
mkdir -p ${gtex_dir}/eqtl/$gtex_eur_eqtl_dir
cd ${gtex_dir}/eqtl/$gtex_eur_eqtl_dir
wget -nc $gtex_eur_eqtl_http
tar -xvf $gtex_eur_eqtl_filename

# download GTEx whole blood TPM
gtex_wb_tpm_https=https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_tpm/gene_tpm_2017-06-05_v8_whole_blood.gct.gz
gtex_wb_tpm=$(basename $gtex_wb_tpm_https)
mkdir -p ${gtex_dir}/tpm_wb
cd ${gtex_dir}/tpm_wb
wget -nc $gtex_wb_tpm_https
gunzip $gtex_wb_tpm
gtex_wb_tpm_name=$(echo $gtex_wb_tpm | sed 's/\.gct.gz$//')
gtex_wb_tpm_name_updtd=${gtex_wb_tpm_name}.tsv
sed '1,2d' ${gtex_wb_tpm_name}.gct > $gtex_wb_tpm_name_updtd
gzip ${gtex_wb_tpm_name}.gct

# download GTEx median TPM per tissue
gtex_median_http=https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
mkdir -p ${gtex_dir}/median_tpm
cd ${gtex_dir}/median_tpm
gtex_median_tpm=$(basename $gtex_median_http)
wget -nc $gtex_median_http
# gunzip and clean file 
gunzip $gtex_median_tpm
gtex_median_name=$(echo $gtex_median_tpm | sed 's/\.gct.gz$//')
gtex_median_name_updtd=${gtex_median_name}.tsv
sed '1,2d' ${gtex_median_name}.gct > $gtex_median_name_updtd
gzip ${gtex_median_name}.gct

### eQTLGen
# File with significant (FDR < 0.05) cis-eQTL results
eqtlgen_https=https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
eqtlgen_dir=${out_dir}/eqtlgen
mkdir -p $eqtlgen_dir
cd $eqtlgen_dir
wget -nc $eqtlgen_https

### Roadmap Epignomics Project ChromHMM 
chromhmm_pbmcs_https=https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E062_15_coreMarks_hg38lift_mnemonics.bed.gz
chromhmm_dir=${out_dir}/chromhmm
mkdir -p $chromhmm_dir
cd $chromhmm_dir
wget -nc $chromhmm_pbmcs_https
chromhmm_pbmcs_filename=$(basename $chromhmm_pbmcs_https)

### xCell 
xcell_https=https://github.com/dviraran/xCell/raw/master/data/xCell.data.rda 
xcell_dir=${out_dir}/xcell
mkdir -p $xcell_dir
cd $xcell_dir
wget -nc $xcell_https 

### GnomAD
gnomad_dir=${out_dir}/gnomad
mkdir -p $gnomad_dir
cd $gnomad_dir
# gnomAD LoF metrics, downloaded from here https://gnomad.broadinstitute.org/downloads/#v2-constraint
gnomad_lof_https=https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz
wget -nc $gnomad_lof_https
gnomad_lof_filename=$(basename $gnomad_lof_https)
gnomad_lof=$(echo $gnomad_lof_filename | sed 's/\.bgz$//')
gnomad_lof_gzip=${gnomad_lof}.gz
cp $gnomad_lof_filename $gnomad_lof_gzip
gunzip $gnomad_lof_gzip
rm $gnomad_lof_gzip

# gnomAD whole genome constraint z-scores
# gnomAD files downloaded from here: https://gnomad.broadinstitute.org/downloads#v3-genomic-constraint
gnomad_zscore_https=https://storage.googleapis.com/gnomad-nc-constraint-v31-paper/download_files/constraint_z_genome_1kb.qc.download.txt.gz
wget -nc $gnomad_zscore_https
gnomad_zscore_file=$(basename $gnomad_zscore_https)
gnomad_zscore_gzip=$(echo $gnomad_zscore_file | sed 's/\.gz$//')
gunzip ${gnomad_zscore_file}
tail -n +2 "$gnomad_zscore_gzip" > ${gnomad_zscore_gzip}.clean

# gnomAD gene sets from Karczewski et al., Nature 2020
mkdir supplement_2020
cd supplement_2020
gnomad_supplement_2020=https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2308-7/MediaObjects/41586_2020_2308_MOESM4_ESM.zip
wget -nc $gnomad_supplement_2020
gnomad_genesets_filename=$(basename $gnomad_supplement_2020)
unzip $gnomad_genesets_filename

### pHaplo and pTriplo from Collins et al., Cell 2022 Supplementary Table 7
phaplo_ptriplo_dir=${out_dir}/phaplo_ptriplo
mkdir -p $phaplo_ptriplo_dir
cd $phaplo_ptriplo_dir
collins_https=https://ars.els-cdn.com/content/image/1-s2.0-S0092867422007887-mmc7.xlsx
wget -nc $collins_https

### EDS from Wang et al., AJHG 2020 Supplementary Table 1 
eds_dir=${out_dir}/eds
mkdir -p $eds_dir
cd $eds_dir
eds_https=https://ars.els-cdn.com/content/image/1-s2.0-S0002929720300124-mmc2.xlsx
wget -nc $eds_https

### Episcore from Han et al., Nature Communications 2018 Supplementary Data 3 
episcore_dir=${out_dir}/episcore
mkdir -p $episcore_dir
cd $episcore_dir
episcore_https=https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-04552-7/MediaObjects/41467_2018_4552_MOESM5_ESM.xlsx
wget -nc $episcore_https

### GERP ++ elements conservation 
gerp_dir=${out_dir}/conservation/gerp
mkdir -p $gerp_dir
cd $gerp_dir
gerp_https=https://bds.mpi-cbg.de/hillerlab/120MammalAlignment/Human120way/data/conservation/gerpElements_hg38_multiz120Mammals.bed.gz
wget -nc $gerp_https
gerp_gzip=$(basename $gerp_https)
gerp_filename=$(echo $gerp_gzip | sed 's/\.gz$//')
gunzip -c $gerp_gzip > $gerp_filename

### PhyloP 100way 
phlop_dir=${out_dir}/conservation/phylop
mkdir -p $phlop_dir
cd $phlop_dir
phylop_https=http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
wget -nc $phylop_https

### Human Accelerated Regions (HARs)
hars_dir=${out_dir}/conservation/hars
mkdir -p $hars_dir
cd $hars_dir
# HARs from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180714
hars_dir_https=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180714/suppl/GSE180714%5FHARs%2Ebed%2Egz
wget -nc $hars_dir_https
zcat GSE180714_HARs.bed.gz | awk 'NR > 1 {OFS="\t"; print $1, $2, $3, $4}' > GSE180714_HARs.bed

### DECIPHER developmental disorder genes 
decipher_dir=${out_dir}/decipher
mkdir -p $decipher_dir
cd $decipher_dir
decipher_https=https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz
wget -nc $decipher_https

### OpenTargets approved drug targets 
opentargets_dir=${out_dir}/opentargets
mkdir -p $opentargets_dir
cd $opentargets_dir
opentargets_https=http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/targets/ 
wget -nc --no-parent --recursive $opentargets_https

### CpG islands
cpg_island_dir=${out_dir}/cpg_islands/
mkdir -p $cpg_island_dir
cd $cpg_island_dir
cpg_island_https=http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz
wget -nc $cpg_island_https
zcat cpgIslandExt.txt.gz | awk 'NR > 1 {OFS="\t"; print $2, $3, $4}' > cpgIslandExt.bed

### ENCODE cCREs 
encode_ccre_dir=${out_dir}/encode/encode_c_cres
mkdir -p $encode_ccre_dir
cd $encode_ccre_dir
# all CTCF-only binding sites
all_ctcf_only_https=https://downloads.wenglab.org/Registry-V3/GRCh38-cCREs.CTCF-only.bed
wget -nc $all_ctcf_only_https
# CD14 monocyte cCREs 
cd14_monocyte_https=https://downloads.wenglab.org/Registry-V3/Seven-Group/ENCFF389PZY_ENCFF587XGD_ENCFF184NWF_ENCFF496PSJ.7group.bed
wget -nc $cd14_monocyte_https
# B-cell cCREs
bcell_ctcf_https=https://downloads.wenglab.org/Registry-V3/Seven-Group/ENCFF035DJL.7group.bed
wget -nc $bcell_ctcf_https
# Neutrophil cCREs
neutrophil_https=https://downloads.wenglab.org/Registry-V3/Seven-Group/ENCFF685DZI_ENCFF311TAY_ENCFF300LXQ.7group.bed
wget -nc $neutrophil_https

### polyA sites database 
polya_site_dir=${out_dir}/polyA_site
mkdir -p $polya_site_dir
cd $polya_site_dir
polya_bed_https=https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz
wget -nc $polya_bed_https
polya_atlas_https=https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.tsv.gz
wget -nc $polya_atlas_https

### 4D nucleome data (downloaded manually) 
# 4DNFILYQ1PAY - A/B compartments, in situ Hi-C on gm12878 with MboI and bio-dATP https://data.4dnucleome.org/files-processed/4DNFILYQ1PAY/
# 4DNFIVK5JOFU - TAD boundaries, in situ Hi-C on gm12878 with MboI and bio-dATP https://data.4dnucleome.org/files-processed/4DNFIVK5JOFU/
# 4DNFIJL18YS3 - TAD boundaries HMEC, in situ Hi-C on CC2551 with MboI and bio-dATP https://data.4dnucleome.org/files-processed/4DNFIJL18YS3/
# 4DNFICLU9GUP - TAD boundaries NHEK, in situ Hi-C on 192627 with MboI and bio-dATP https://data.4dnucleome.org/files-processed/4DNFICLU9GUP/
# 4DNFI9MZWZF7 - TAD boundaries HUVEC, in situ Hi-C on CC2517 with MboI and bio-dATP https://data.4dnucleome.org/files-processed/4DNFI9MZWZF7/
# 4DNFIMNT2VYL - TAD boundaries IMR90, in situ Hi-C on IMR90 with MboI and bio-dATP https://data.4dnucleome.org/files-processed/4DNFIMNT2VYL/

### OMIM - data request required (downloaded manually)

### COSMIC login required - downloaded from here: https://cancer.sanger.ac.uk/census
