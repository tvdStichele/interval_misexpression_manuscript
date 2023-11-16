#!/bin/python
import pandas as pd 
import pysam 
from io import StringIO
from pybedtools import BedTool
from pathlib import Path
import argparse

def main():
    # arguments passed 
    parser = argparse.ArgumentParser(description="")
    # required
    parser.add_argument("--gene_id", 
                        help ="Ensembl ID of gene", 
                        required = True)
    parser.add_argument("--chrom", 
                        help ="chromosome", 
                        required = True)
    parser.add_argument("--vcf", 
                        help = "VCF file", 
                        required = True)
    parser.add_argument("--rna_ids", 
                        help="List of RNA IDs passing QC", 
                        required=True)
    parser.add_argument("--paired_smpls", 
                        help = "file with paired EGAN and RNA sample IDs", 
                        required = True)
    parser.add_argument("--intersect_bed", 
                        help = "Bed file with intersecting genes and variants", 
                        required = True)
    parser.add_argument("--root", 
                        help = "root directory output", 
                        required = True)
    # set arguments to global variable 
    global args
    args = parser.parse_args()
    # load arguments 
    gene_id = args.gene_id
    chrom = args.chrom
    vcf_path = args.vcf
    smpl_id_pass_qc_set = set(pd.read_csv(args.rna_ids, sep=",", header=None)[0])
    wgs_rna_paired_smpls_path = args.paired_smpls
    intersect_bed_path = args.intersect_bed
    root_dir = args.root 
    
    root_dir_path = Path(root_dir)
    root_dir_path.mkdir(parents=True, exist_ok=True)
    # load intersection bed file 
    intersect_bed_df = pd.read_csv(intersect_bed_path, sep="\t")
    # extract entries with gene
    gene_id_intersect_bed_df = intersect_bed_df[intersect_bed_df.gene_id == gene_id]
    # generate set of variant IDs inside gene windows
    gene_id_intersect_no_dupl_df = gene_id_intersect_bed_df[["chrom_vrnt", "vrnt_id", "start_vrnt", "end_vrnt"]].drop_duplicates()
    intersect_vrnt_ids_list = gene_id_intersect_no_dupl_df.vrnt_id.unique()
    print(f"{chrom} number of variants in {gene_id} windows: {len(intersect_vrnt_ids_list)}")

    vcf = pysam.VariantFile(vcf_path, mode = "r")
    # load RNA-EGAN ID linker and subset to samples passing QC
    wgs_rna_paired_smpls_df = pd.read_csv(wgs_rna_paired_smpls_path, sep="\t")
    egan_ids_with_rna = wgs_rna_paired_smpls_df[wgs_rna_paired_smpls_df.rna_id.isin(smpl_id_pass_qc_set)].egan_id.tolist()
    # subset vcf to samples with expression data 
    vcf_samples = [sample for sample in vcf.header.samples]
    vcf_egan_ids_with_rna = set(egan_ids_with_rna).intersection(set(vcf_samples))
    vcf.subset_samples(vcf_egan_ids_with_rna)
    vcf_samples_with_rna = [sample for sample in vcf.header.samples]
    print(f"Number of samples in VCF with RNA ID and passing QC: {len(vcf_samples_with_rna)}")
    # subset egan ID and RNA ID linker file to samples with genotype calls and passing QC 
    wgs_rna_paired_smpls_with_gts_df = wgs_rna_paired_smpls_df[wgs_rna_paired_smpls_df.egan_id.isin(vcf_samples_with_rna)]
    rna_id_pass_qc_with_gts = wgs_rna_paired_smpls_with_gts_df.rna_id.unique().tolist()
    print(f"Number of RNA IDs passing QC with paired sample in VCF: {len(rna_id_pass_qc_with_gts)}")

    vrnts_gt_intersect_dir = root_dir_path.joinpath(f"vrnts_gts_intersect/{chrom}/")
    Path(vrnts_gt_intersect_dir).mkdir(parents=True, exist_ok=True)
    vrnts_annotate_path = vrnts_gt_intersect_dir.joinpath(f"{chrom}_{gene_id}_vrnts_gts_intersect.tsv")

    all_gts_set = {(1, 0), (0, 1), (1, 1), (0, 0), (None, None)}
    carriers_lt_50perc = {(1, 0), (0, 1), (1, 1)}
    carriers_gt_50perc = {(1, 0), (0, 1), (0, 0)}
    with open(vrnts_annotate_path, "w") as f_out: 
        header='\t'.join(['vrnt_id', 'egan_id', 'genotype', 'AF', 'vrnt_type'])
        f_out.write(f"{header}\n")
        for vrnt_id in intersect_vrnt_ids_list:
            # get chromosome, start and end of variant
            chrom_vrnt, pos, end = [gene_id_intersect_no_dupl_df[gene_id_intersect_no_dupl_df.vrnt_id == vrnt_id][col].item() for col in ["chrom_vrnt", "start_vrnt", "end_vrnt"]]
            records = vcf.fetch(f"chr{chrom_vrnt}", pos, end)
            found_vrnt_id = False
            for rec in records: 
                chrom, pos, ref, alt_alleles = rec.chrom, rec.pos, rec.ref, rec.alts
                # check for multiallelic alleles
                if len(alt_alleles) > 1: 
                    raise ValueError(f"Multiallelic entry in vcf at: {chrom}, {pos}")
                else: 
                    alt = alt_alleles[0]
                vcf_vrnt_id = f"{chrom}:{pos}:{ref}:{alt}"
                if vrnt_id == vcf_vrnt_id:
                    found_vrnt_id = True 
                    gts = [s["GT"] for s in rec.samples.values()]
                    gts_set = set(gts)
                    # check genotypes 
                    if not gts_set <= all_gts_set: 
                        raise ValueError(f"Encountered different genotype at {chrom}, {pos}: {gts_set}.")
                    # assign carrier genotypes 
                    af = rec.info["AF"]
                    if af < 0.5: 
                        carrier_gts = carriers_lt_50perc  
                    else: 
                        carrier_gts = carriers_gt_50perc  
                    # write if variant has carriers 
                    if len(carrier_gts.intersection(gts_set)) > 0:
                        if len(ref) == 1 and len(alt) == 1: 
                            vrnt_type = "snp"
                        else:
                            vrnt_type = "indel"
                        gts_dict = {i:gt for i, gt in enumerate(gts) if gt in carrier_gts} 
                        for i in gts_dict.keys():
                            line= '\t'.join([vrnt_id, vcf.header.samples[i], str(gts_dict[i]), str(af), vrnt_type])
                            f_out.write(f"{line}\n")
            if not found_vrnt_id: 
                raise ValueError(f"Did not find {vrnt_id} in {vcf_path}")
    f_out.close()
    vcf.close()
    
    
if __name__ == "__main__":
    main()
    