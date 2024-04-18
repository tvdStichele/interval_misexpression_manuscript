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
    parser = argparse.ArgumentParser(description="")
    # required 
    parser.add_argument("--express", 
                        help="Flat gene expression matrix. CSV format. Columns:", 
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
                        help = "File with paired RNA and WGS sample IDs. RNA IDs in first column. WGS IDs in second.", 
                        required = True)
    parser.add_argument("--vep_msc", 
                        help = "File with variant most severe consequences", 
                        required = True)
    parser.add_argument("--window", 
                        help = "Size of window around gene in base pairs", 
                        type=int,
                        required = True)
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
    vep_msc_path = args.vep_msc
    window = args.window
    root_dir = args.root
    
    print("Inputs:")
    print(f"- Chromosome: {chrom}")
    print(f"- Window size: {window}bp")
    print("")
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
                gene_id_list, chrom_list, start_list, end_list = [], [], [], []
                gene_id_list.append(gtf.gene_id.split('.')[0])
                chrom_list.append(gtf.contig)
                start_list.append(gtf.start)
                end_list.append((gtf.end))
                # check or write output
                if len(gene_id_list) > 1 or len(chrom_list) > 1 or len(start_list) > 1 or len(end_list) > 1:
                    print(f"{gene_id} has multiple entries in gencode file - excluded from output file.")
                elif len(chrom_list) == 0 or len(start_list) == 0 or len(end_list) == 0: 
                    print(f"{gene_id} has no entries in gencode file - excluded from output file.")
                else: 
                    gene_id, gtf_chrom, start, end = gene_id_list[0], chrom_list[0], start_list[0], end_list[0]
                    gene_id_pass_qc_on_chrom.append(gene_id)
                    chrom_num = gtf_chrom.split("chr")[1]
                    f.write(f"{chrom_num}\t{start}\t{end}\t{gene_id}\n")
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
    
    ### intersect variants with gene windows 
    intersect_bed_dir = root_dir_path.joinpath("intersect_bed")
    intersect_bed_dir.mkdir(parents=True, exist_ok=True)
    intersect_bed_path = intersect_bed_dir.joinpath(f"{chrom}_intersect.bed")

    vrnts_bed = BedTool(vrnts_bed_path)
    gencode_bed= BedTool(gene_bed_path)

    intersect_bed_str = StringIO(str(gencode_bed.window(vrnts_bed, w=window)))
    columns={0:"chrom_gene", 1:"start_gene", 2:"end_gene", 3: "gene_id", 
             4:"chrom_vrnt", 5:"start_vrnt", 6:"end_vrnt", 7:"vrnt_id"}
    intersect_bed_df = pd.read_csv(intersect_bed_str, 
                                   sep="\t", 
                                   header=None).rename(columns=columns)
    intersect_bed_df.to_csv(intersect_bed_path, sep="\t", index=False, header=None)
    # generate set of variant IDs inside gene windows
    intersect_vrnt_ids_no_dupl_df = intersect_bed_df[["chrom_vrnt", "vrnt_id", "start_vrnt", "end_vrnt"]].drop_duplicates()
    intersect_vrnt_ids_list = intersect_vrnt_ids_no_dupl_df.vrnt_id.unique()
    num_intersect_vrnt_ids = len(intersect_vrnt_ids_list)
    print(f"Number of variants intersecting gene windows: {num_intersect_vrnt_ids}")
    
    ### Genotypes for all variant IDs in windows
    count = 0
    vrnt_gt_egan_dict = {}
    for vrnt_id in intersect_vrnt_ids_list:
        # get chromosome, start and end of SV 
        chrom_vrnt, pos, end, = [intersect_vrnt_ids_no_dupl_df[intersect_vrnt_ids_no_dupl_df.vrnt_id == vrnt_id][col].item() for col in ["chrom_vrnt", "start_vrnt", "end_vrnt"]]
        # search VCF for variant ID
        records = vcf.fetch(str(chrom_vrnt), pos-1, end)
        found_vrnt_id = False
        for rec in records: 
            vcf_vrnt_id = str(rec.id)
            if vrnt_id == vcf_vrnt_id:
                found_vrnt_id = True 
                # collect genotypes
                gts = [s["GT"] for s in rec.samples.values()]
                for i, gt in enumerate(gts): 
                    vrnt_gt_egan_dict[count] = [vrnt_id, vcf_samples_with_rna[i], gt]
                    count += 1 
        if not found_vrnt_id: 
            raise ValueError(f"Did not find {vrnt_id} in {vcf_path}")
    vrnt_gt_egan_nogene_df = pd.DataFrame.from_dict(vrnt_gt_egan_dict, 
                                                    orient="index", 
                                                    columns=["vrnt_id", "egan_id", "genotype"])
    if vrnt_gt_egan_nogene_df.shape[0] != num_intersect_vrnt_ids * len(rna_id_pass_qc_sv_calls):
        raise ValueError("Number of variant genotypes does not much variant number by samples.")
    
    ### add intersecting gene, SV info, VEP MSC and expression
    # load SV info
    sv_info_df = pd.read_csv(sv_info_path, sep="\t", dtype={"plinkID":str})
    sv_info_id_af_df = sv_info_df[["plinkID", "AF", "SVTYPE"]].rename(columns={"plinkID":"vrnt_id"})
    # load most severe consequence 
    vep_msc_df = pd.read_csv(vep_msc_path, sep="\t", dtype={"vrnt_id": str})
    msc_list = vep_msc_df.Consequence.unique().tolist()
    vep_msc_cnsqn_df = vep_msc_df.rename(columns={"Uploaded_variation": "vrnt_id"})[["vrnt_id", "Consequence"]]
    # merge 
    dfs_to_merge = [vrnt_gt_egan_nogene_df, 
                    intersect_bed_df[["vrnt_id","gene_id"]].drop_duplicates(),
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
    # add window name
    sv_intersect_express_info_df["window"] = f"gene_body_{window}"
    # only keep carriers 
    sv_intersect_express_info_carrier_df = sv_intersect_express_info_df[sv_intersect_express_info_df.genotype.isin([(0, 1), (1, 1)])]
    # write genotypes for SVs  
    express_gt_info_dir = root_dir_path.joinpath("express_carrier_info")
    express_gt_info_dir.mkdir(parents=True, exist_ok=True)
    intersect_bed_path = express_gt_info_dir.joinpath(f"{chrom}_express_carrier_info.tsv")
    sv_intersect_express_info_carrier_df.to_csv(intersect_bed_path, sep="\t", index=False)
    
    
if __name__ == "__main__":
    main()