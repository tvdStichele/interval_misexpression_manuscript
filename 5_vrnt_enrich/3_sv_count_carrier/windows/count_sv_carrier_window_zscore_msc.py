#!/bin/python
import pandas as pd 
import numpy as np
import pysam 
import argparse
from pybedtools import BedTool
from pathlib import Path
from io import StringIO
from functools import reduce

def main():
    # arguments passed 
    parser = argparse.ArgumentParser(description="Variants correlated with misexpression events")
    # required 
    parser.add_argument("--express", 
                        help="Flat gene expression matrix.", 
                        required=True)
    parser.add_argument("--chrom", 
                        help ="chromosome", 
                        required = True)
    parser.add_argument("--sv_info", 
                        help = "structural variants information file", 
                        required = True)
    parser.add_argument("--vcf", 
                        help = "VCF file", 
                        required = True)
    parser.add_argument("--gencode", 
                        help = "gencode file", 
                        required = True)
    parser.add_argument("--paired_smpls", 
                        help = "file with paired EGAN and RNA sample IDs", 
                        required = True)
    parser.add_argument("--vep_msc", 
                        help = "File with variant most severe consequences", 
                        required = True)
    parser.add_argument('--z_cutoff', 
                        nargs='+', 
                        type=float, 
                        required = True,
                        help = "list of TPM cutoffs")
    parser.add_argument("--window_start", 
                        help = "Start position of upstream and downstream window relative to gene start and end", 
                        type=int,
                        required = True)
    parser.add_argument("--window_step_size", 
                        help = "Window step-size", 
                        type=int,
                        required = True)
    parser.add_argument("--window_max", 
                        help = "Window step-size", 
                        type=int,
                        required = True)
    parser.add_argument("--af_range", 
                        nargs='+', 
                        type=float,
                        required = True, 
                        help = "Allele frequency range 0 to 1")
    parser.add_argument("--root", 
                        help = "root directory output", 
                        required = True)
    # set arguments to global variable 
    global args
    args = parser.parse_args()
    # load arguments 
    ge_matrix_flat_path = args.express
    chrom = args.chrom
    sv_info_path = args.sv_info
    vcf_path=args.vcf
    gencode_path=args.gencode
    wgs_rna_paired_smpls_path = args.paired_smpls
    z_cutoff_list = args.z_cutoff
    window_start = args.window_start
    window_step_size = args.window_step_size
    window_max = args.window_max
    af_lower, af_upper = args.af_range
    vep_msc_path = args.vep_msc
    root_dir = args.root

    print("Inputs:")
    print(f"- Chromosome: {chrom}")
    print(f"- AF range: {af_lower}-{af_upper}")
    print(f"- Window start: {window_start}")
    print(f"- Window size: {window_step_size}")
    print(f"- Max window: {window_max}")
    print(f"- Z-score cutoffs: {z_cutoff_list}")
    print("")
    
    # constants 
    tpm_cutoff = 0.5 
    
    # create root directory 
    root_dir_path = Path(root_dir)
    root_dir_path.mkdir(parents=True, exist_ok=True)
    
    ### Collect samples with VCF calls and RNA sequencing 
    # read in gene expression file 
    ge_matrix_flat_df = pd.read_csv(ge_matrix_flat_path, sep=",")
    gene_id_pass_qc_set = set(ge_matrix_flat_df.gene_id.unique())
    print(f"Gene IDs passing filters: {len(gene_id_pass_qc_set)}")
    smpl_id_pass_qc_set = set(ge_matrix_flat_df.rna_id.unique())
    print(f"RNA-seq sample IDs passing QC: {len(smpl_id_pass_qc_set)}")
    print("")
    # egan ID, RNA ID sample links 
    wgs_rna_paired_smpls_df = pd.read_csv(wgs_rna_paired_smpls_path, sep="\t")
    egan_ids_with_rna = wgs_rna_paired_smpls_df[wgs_rna_paired_smpls_df.rna_id.isin(smpl_id_pass_qc_set)].egan_id.tolist()
    # load VCF and subset to EGAN IDs with RNA 
    print("Loading input VCF ...")
    vcf_path = vcf_path
    vcf = pysam.VariantFile(vcf_path, mode = "r")
    print("VCF loaded.")
    print("Subset VCF to samples with RNA-seq ...")
    vcf_samples = [sample for sample in vcf.header.samples]
    vcf_egan_ids_with_rna = set(egan_ids_with_rna).intersection(set(vcf_samples))
    vcf.subset_samples(vcf_egan_ids_with_rna)
    vcf_samples_with_rna = [sample for sample in vcf.header.samples]
    print(f"Number of samples in VCF with RNA ID and passing QC: {len(vcf_samples_with_rna)}")
    # subset egan ID and RNA ID links to samples with SV calls and passing QC 
    wgs_rna_paired_smpls_with_sv_calls_df = wgs_rna_paired_smpls_df[wgs_rna_paired_smpls_df.egan_id.isin(vcf_samples_with_rna)]
    # write EGAN-RNA ID pairs to file
    egan_rna_smpls_dir = root_dir_path.joinpath("egan_rna_smpls")
    Path(egan_rna_smpls_dir).mkdir(parents=True, exist_ok=True)
    wgs_rna_paired_smpls_with_sv_calls_df.to_csv(egan_rna_smpls_dir.joinpath("egan_rna_ids_paired_pass_qc.tsv"), sep="\t", index=False)
    rna_id_pass_qc_sv_calls = wgs_rna_paired_smpls_with_sv_calls_df.rna_id.unique().tolist()
    print(f"Number of RNA IDs passing QC: {len(rna_id_pass_qc_sv_calls)}")
    
    ### write bed file for genes on chromosome passing QC 
    gene_bed_dir = root_dir_path.joinpath("genes_bed")
    gene_bed_dir.mkdir(parents=True, exist_ok=True)
    gene_bed_path = gene_bed_dir.joinpath(f"{chrom}_genes.bed")
    gene_id_pass_qc_on_chrom = []
    with open(gene_bed_path, "w") as f:
        for gtf in pysam.TabixFile(gencode_path).fetch(chrom, parser = pysam.asGTF()):
            if gtf.feature == "gene" and gtf.gene_id.split('.')[0] in gene_id_pass_qc_set:
                # check for multiple entries with same name 
                gene_id_list, chrom_list, start_list, end_list, strand_list = [], [], [], [], []
                gene_id_list.append(gtf.gene_id.split('.')[0])
                chrom_list.append(gtf.contig)
                start_list.append(gtf.start)
                end_list.append((gtf.end))
                strand_list.append((gtf.strand))
                # check or write output
                if len(gene_id_list) > 1 or len(chrom_list) > 1 or len(start_list) > 1 or len(end_list) > 1:
                    print(f"{gene_id} has multiple entries in gencode file - excluded from output file.")
                elif len(chrom_list) == 0 or len(start_list) == 0 or len(end_list) == 0: 
                    print(f"{gene_id} has no entries in gencode file - excluded from output file.")
                else: 
                    gene_id, gtf_chrom, start, end, strand = gene_id_list[0], chrom_list[0], start_list[0], end_list[0], strand_list[0]
                    gene_id_pass_qc_on_chrom.append(gene_id)
                    chrom_num = gtf_chrom.split("chr")[1]
                    f.write(f"{chrom_num}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n")
    print(f"Number of genes pass QC on {chrom}: {len(gene_id_pass_qc_on_chrom)}")
    
    # subset gene expression file to genes on chromosome 
    ge_matrix_flat_chrom_df = ge_matrix_flat_df[ge_matrix_flat_df.gene_id.isin(gene_id_pass_qc_on_chrom)]
    # subset gene expression file to samples with SV calls 
    ge_matrix_flat_chrom_egan_df = pd.merge(ge_matrix_flat_chrom_df, wgs_rna_paired_smpls_with_sv_calls_df, how="inner", on="rna_id")
    print(f"Sample IDs in gene expression matrix with SV calls: {len(ge_matrix_flat_chrom_egan_df.egan_id.unique())}")
    # write to file 
    ge_matrix_flat_chrom_egan_dir = root_dir_path.joinpath("express_mtx")
    Path(ge_matrix_flat_chrom_egan_dir).mkdir(parents=True, exist_ok=True)
    ge_matrix_flat_chrom_egan_df.to_csv(ge_matrix_flat_chrom_egan_dir.joinpath(f"{chrom}_ge_mtx_flat.tsv"), sep="\t", index=False)
    # add gene-sample pair to expression dataframe
    ge_matrix_flat_chrom_egan_df["gene_smpl_pair"] = ge_matrix_flat_chrom_egan_df.gene_id + "," + ge_matrix_flat_chrom_egan_df.rna_id
    
    ### write SVs on chromosome to bed file 
    vrnts_bed_dir = root_dir_path.joinpath("vrnts_bed")
    vrnts_bed_dir.mkdir(parents=True, exist_ok=True)
    vrnts_bed_path = vrnts_bed_dir.joinpath(f"{chrom}_vrnts.bed")
    with open(sv_info_path, "r") as f_in, open(vrnts_bed_path, "w") as f_out:
        for line in f_in:
            if line.startswith("plinkID"): 
                continue
            else: 
                vrnt_id, sv_chrom, pos, end = line.split("\t")[0], line.split("\t")[2], line.split("\t")[3], line.split("\t")[4]
                if sv_chrom == chrom:
                    chrom_num = sv_chrom.split("chr")[1]
                    f_out.write(f"{chrom_num}\t{pos}\t{end}\t{vrnt_id}\n")
    
    ### SV information 
    sv_info_df = pd.read_csv(sv_info_path, sep="\t", dtype={"plinkID":str})
    sv_types_list = sv_info_df.SVTYPE.unique().tolist()
    print(f"SV types in SV info file: {sv_types_list}")
    sv_info_id_af_df = sv_info_df[["plinkID", "AF", "SVTYPE"]].rename(columns={"plinkID":"vrnt_id"})

    ### generate output directories 
    # variant-window intersection directory 
    intersect_bed_dir=root_dir_path.joinpath("intersect_bed")
    Path(intersect_bed_dir).mkdir(parents=True, exist_ok=True)
    # genotypes of variants in windows directory 
    intersect_vrnt_gts_dir = root_dir_path.joinpath("intersect_vrnt_gts")
    Path(intersect_vrnt_gts_dir).mkdir(parents=True, exist_ok=True)

    # load variant and gene bed files 
    vrnts_bed = BedTool(vrnts_bed_path)
    genes_bed = BedTool(gene_bed_path)
    # column names for variant window intersection bed file 
    intersect_bed_columns={0:"chrom_gene", 1:"start_gene", 2:"end_gene", 3: "gene_id", 4: "score", 5: "strand", 
                           6:"chrom_vrnt", 7:"start_vrnt", 8:"end_vrnt", 9:"vrnt_id"}

    # load most severe consequence 
    vep_msc_df = pd.read_csv(vep_msc_path, sep="\t", dtype={"vrnt_id": str})
    msc_list = vep_msc_df.Consequence.unique().tolist()
    vep_msc_cnsqn_df = vep_msc_df.rename(columns={"Uploaded_variation": "vrnt_id"})[["vrnt_id", "Consequence"]]
    
    entry = 0 
    carrier_count = {}
    for direction in ["upstream", "downstream"]:
        window_size = window_start
        seen_sv_gene_pairs = set()
        while window_size <= window_max:
            print(f"Counting variants {direction}, {window_size}bp ...")
            if direction == "upstream": 
                l, r = window_size, 0 
            else: 
                l, r = 0, window_size
            vrnts_window_intersect_str = StringIO(str(genes_bed.window(vrnts_bed, 
                                                                       l=l, 
                                                                       r=r,
                                                                       sw=True
                                                                      )))        
            intersect_bed_df = pd.read_csv(vrnts_window_intersect_str, sep="\t", header=None).rename(columns=intersect_bed_columns)
            # variant gene ID column 
            intersect_bed_df["vrnt_gene_id"] = intersect_bed_df["vrnt_id"].astype(str) + "," + intersect_bed_df["gene_id"].astype(str)
            # subset to variant gene pairs that have not been seen in previous windows    
            sv_gene_pairs_in_window = set(intersect_bed_df.vrnt_gene_id.unique())
            new_sv_gene_pairs = sv_gene_pairs_in_window - seen_sv_gene_pairs
            if new_sv_gene_pairs.union(seen_sv_gene_pairs) != sv_gene_pairs_in_window: 
                raise ValueError("Missing variants from windows.")
            intersect_bed_new_sv_gene_df = intersect_bed_df[intersect_bed_df["vrnt_gene_id"].isin(new_sv_gene_pairs)] 
            # generate set of variant IDs inside gene windows
            intersect_vrnt_ids_no_dupl_df = intersect_bed_new_sv_gene_df[["chrom_vrnt", "vrnt_id", "start_vrnt", "end_vrnt"]].drop_duplicates()
            intersect_vrnt_ids_list = intersect_vrnt_ids_no_dupl_df.vrnt_id.unique()
            # genotypes for variant IDs in windows
            vrnt_gt_egan_dict = {}
            count = 0
            for vrnt_id in intersect_vrnt_ids_list:
                # get chromosome, start and end of SV 
                chrom_vrnt, pos, end, = [intersect_vrnt_ids_no_dupl_df[intersect_vrnt_ids_no_dupl_df.vrnt_id == vrnt_id][col].item() for col in ["chrom_vrnt", "start_vrnt", "end_vrnt"]]
                records = vcf.fetch(str(chrom_vrnt), pos-1, end)
                found_vrnt_id = False
                for rec in records: 
                    vcf_vrnt_id = str(rec.id)
                    if vrnt_id == vcf_vrnt_id:
                        found_vrnt_id = True 
                        gts = [s["GT"] for s in rec.samples.values()]
                        for i, gt in enumerate(gts): 
                            vrnt_gt_egan_dict[count] = [vrnt_id, vcf.header.samples[i], gt]
                            count += 1 
                if not found_vrnt_id: 
                    raise ValueError(f"Did not find {vrnt_id} in {vcf_path}")
            vrnt_gt_egan_nogene_df = pd.DataFrame.from_dict(vrnt_gt_egan_dict, orient="index", columns=["vrnt_id", "egan_id", "genotype"])
            vrnt_gt_egan_nogene_df = vrnt_gt_egan_nogene_df.astype({"genotype":str})

            # merge 
            dfs_to_merge = [vrnt_gt_egan_nogene_df, 
                            intersect_bed_new_sv_gene_df[["vrnt_id", "gene_id","vrnt_gene_id"]].drop_duplicates(),
                            sv_info_id_af_df, 
                            vep_msc_cnsqn_df
                           ]
            df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['vrnt_id'],
                                                        how='inner'), dfs_to_merge)
            # add expression 
            sv_intersect_express_info_df = pd.merge(df_merged, 
                                                 ge_matrix_flat_chrom_egan_df,                  
                                                 how="inner",
                                                 on=["egan_id", "gene_id"])
            # write to file 
            intersect_vrnt_gts_path = intersect_vrnt_gts_dir.joinpath(f"{chrom}_intersect_vrnt_gts_{direction}_{window_size}.tsv")
            sv_intersect_express_info_df.to_csv(intersect_vrnt_gts_path, sep="\t", index=False)

            ### count carriers in test and controls across different Z-score cutoffs 
            for z_cutoff in z_cutoff_list: 
                carrier_count[entry] = [chrom, direction, window_size, z_cutoff]
                ### get control and test gene-sample pairs based on TPM cutoff 

                # get gene-sample pairs passing z-score and TPM cutoff 
                misexp_df = ge_matrix_flat_chrom_egan_df[(ge_matrix_flat_chrom_egan_df["z-score"] > z_cutoff) &
                                                         (ge_matrix_flat_chrom_egan_df["TPM"] > tpm_cutoff)
                                                        ]
                test_gene_ids = misexp_df.gene_id.unique()
                test_gene_smpl_pairs = misexp_df.gene_smpl_pair.unique()
                # get gene-sample pairs in control group - same gene set, samples not misexpressing
                cntrl_gene_smpl_pairs = ge_matrix_flat_chrom_egan_df[~ge_matrix_flat_chrom_egan_df.gene_smpl_pair.isin(test_gene_smpl_pairs) &
                                                                     ge_matrix_flat_chrom_egan_df.gene_id.isin(test_gene_ids)
                                                                    ].gene_smpl_pair.unique()
                # check overlap between sets is empty 
                if len(set(test_gene_smpl_pairs).intersection(set(cntrl_gene_smpl_pairs))) != 0:
                    raise ValueError("Overlap between control and test gene-sample pair sets")
                # check number of test and control gene-pairs matches expected 
                total_test_control_pairs = len(test_gene_smpl_pairs) + len(cntrl_gene_smpl_pairs)
                test_genes_by_smpls = len(test_gene_ids) * len(rna_id_pass_qc_sv_calls)
                if total_test_control_pairs !=  test_genes_by_smpls: 
                    raise ValueError("Number of test and control gene pairs does not match test gene IDs by samples")
                carrier_count[entry] += [len(test_gene_ids), len(rna_id_pass_qc_sv_calls), len(test_gene_smpl_pairs), len(cntrl_gene_smpl_pairs)]

                ### count carriers in control and test groups 
                # create column with gene-sample pairs  
                sv_intersect_express_info_df["gene_smpls_id"] = sv_intersect_express_info_df.gene_id.astype(str) + "," + sv_intersect_express_info_df.rna_id.astype(str)
                # get misexpression carriers 
                misexpress_carriers_df = sv_intersect_express_info_df[(sv_intersect_express_info_df.gene_smpls_id.isin(test_gene_smpl_pairs)) & 
                                                                      (sv_intersect_express_info_df.AF >= af_lower) & 
                                                                      (sv_intersect_express_info_df.AF < af_upper) &
                                                                      (sv_intersect_express_info_df.genotype.isin(['(0, 1)', '(1, 1)']))
                                                                     ].copy()
                # get non-misexpression carriers
                cntrl_carriers_df = sv_intersect_express_info_df[(sv_intersect_express_info_df.gene_smpls_id.isin(cntrl_gene_smpl_pairs)) &
                                                                          (sv_intersect_express_info_df.AF >= af_lower) & 
                                                                          (sv_intersect_express_info_df.AF < af_upper) &
                                                                          (sv_intersect_express_info_df.genotype.isin(['(0, 1)', '(1, 1)']))
                                                                         ].copy()
                # add number of carriers 
                gene_smpl_misexpress_carriers = misexpress_carriers_df.gene_smpls_id.unique()
                gene_smpl_cntrl_carriers = cntrl_carriers_df.gene_smpls_id.unique()
                carrier_count[entry] += [len(gene_smpl_misexpress_carriers), len(gene_smpl_cntrl_carriers)]
                # check overlap between sets is empty 
                if len(set(gene_smpl_misexpress_carriers).intersection(set(gene_smpl_cntrl_carriers))) != 0:
                    raise ValueError("Overlap between control and test gene-sample pair carrier sets")

                ### count carriers for SV types and most severe consequence 
                # question remains whether to assign priority to different SV classes
                for sv_type in sv_types_list:
                    misexpress_carriers_sv_type_df = misexpress_carriers_df[misexpress_carriers_df.SVTYPE == sv_type]
                    cntrl_carriers_sv_type_df = cntrl_carriers_df[cntrl_carriers_df.SVTYPE == sv_type]
                    sv_type_misexp_carriers = len(misexpress_carriers_sv_type_df.gene_smpls_id.unique())
                    sv_type_cntrl_carriers = len(cntrl_carriers_sv_type_df.gene_smpls_id.unique())
                    carrier_count[entry] += [sv_type_misexp_carriers, sv_type_cntrl_carriers]
                    for msc in msc_list: 
                        msc_type_misexp_carriers = len(misexpress_carriers_sv_type_df[misexpress_carriers_sv_type_df.Consequence == msc].gene_smpls_id.unique())
                        msc_type_cntrl_carriers = len(cntrl_carriers_sv_type_df[cntrl_carriers_sv_type_df.Consequence == msc].gene_smpls_id.unique())
                        carrier_count[entry] += [msc_type_misexp_carriers, msc_type_cntrl_carriers]
                entry += 1
            # add variant gene pairs to seen set 
            seen_sv_gene_pairs = seen_sv_gene_pairs.union(sv_gene_pairs_in_window)
            window_size += window_step_size

    # write results 
    vrnt_type_cols = []
    for sv_type in sv_types_list: 
        vrnt_type_cols += [f"{sv_type}_misexp", f"{sv_type}_contrl"]
        for msc in msc_list:
            vrnt_type_cols += [f"{sv_type}_{msc}_misexp", f"{sv_type}_{msc}_contrl"]

    columns = ["chrom", "direction", "window_size", "z_cutoff", "misexp_genes", "smpls_pass_qc", "total_misexp", "total_control", "all_sv_misexp", "all_sv_contrl"] + vrnt_type_cols
    carrier_count_df = pd.DataFrame.from_dict(carrier_count, orient="index", columns=columns)
    carrier_count_dir = f"{root_dir}/carrier_count"
    Path(carrier_count_dir).mkdir(parents=True, exist_ok=True)
    carrier_count_path = f"{carrier_count_dir}/{chrom}_carrier_count.tsv"
    carrier_count_df.to_csv(carrier_count_path, sep="\t", index=False)


    
if __name__ == "__main__":
    main()
    
    
    