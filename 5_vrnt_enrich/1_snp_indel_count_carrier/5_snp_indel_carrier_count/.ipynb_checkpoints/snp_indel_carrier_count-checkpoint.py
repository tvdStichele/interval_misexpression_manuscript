import pysam
import glob
import argparse
import pandas as pd 
from pybedtools import BedTool
from pathlib import Path
from io import StringIO

def main():
    # arguments passed 
    parser = argparse.ArgumentParser(description="")
    # required 
    parser.add_argument("--express", 
                        help="flat gene expression matrix.", 
                        required=True)
    parser.add_argument("--chrom", 
                        help ="chromosome", 
                        required = True)
    parser.add_argument("--vcf", 
                        help = "VCF file", 
                        required = True)
    parser.add_argument("--gts_root", 
                        help = "genotypes root", 
                        required = True)
    parser.add_argument("--paired_smpls", 
                        help = "Paired samples", 
                        required = True)
    parser.add_argument("--genes_bed", 
                        help = "Gene bed file", 
                        required = True)
    parser.add_argument("--vrnts_bed", 
                        help = "Variants bed file", 
                        required = True)
    parser.add_argument("--out", 
                        help = "output directory", 
                        required = True)
    # set arguments to global variable 
    global args
    args = parser.parse_args()
    # load input variables
    ge_mtx_path = args.express
    chrom = args.chrom
    vcf_path = args.vcf
    gt_info_root = args.gts_root
    wgs_rna_paired_smpls_path = args.paired_smpls
    genes_bed_path = args.genes_bed
    vrnts_bed_path = args.vrnts_bed
    out_dir = args.out
    
    # constants 
    z_cutoff_list = [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    window_start = 0
    gene_window_size = 10000
    window_step_size = 200000
    window_max = 1000000
    maf_range_list = [[0, 0.01], [0.01, 0.05], [0.05, 0.1], [0.1, 0.25]]
    tpm_cutoff = 0.5
   
    # create output directory
    out_dir_path = Path(out_dir)
    out_dir_path.mkdir(parents=True, exist_ok=True)
    # genotype information directory
    gt_info_root_path = Path(gt_info_root)
    # load gene expression matrix 
    ge_mtx_df = pd.read_csv(ge_mtx_path)

    ### subset gene expresstion matrix to genes on chromsome 
    misexp_genes = pd.read_csv(genes_bed_path, sep="\t", header=None)[3].unique()
    num_misexp_genes = len(misexp_genes)
    print(f"Number of misexpressed genes on chromosome: {num_misexp_genes}")
    #gene_chrom_df = gencode_key_info_df[["gene_id_no_vrsn", "chrom"]].rename(columns={"gene_id_no_vrsn": "gene_id"})
    #ge_mtx_add_chrom_df = pd.merge(ge_mtx_df, gene_chrom_df, on="gene_id", how="inner")
    ge_mtx_chrom_df = ge_mtx_df[ge_mtx_df.gene_id.isin(misexp_genes)]

    ### subset to gene expression matrix to samples with genotyping data 
    # get EGAN IDs with RNA data 
    wgs_rna_paired_smpls_df = pd.read_csv(wgs_rna_paired_smpls_path, sep="\t")
    egan_ids_with_rna = wgs_rna_paired_smpls_df.egan_id.unique()
    # load vcf
    vcf = pysam.VariantFile(vcf_path, mode = "r")
    vcf_samples = [sample for sample in vcf.header.samples]
    vcf_egan_ids_with_rna = set(egan_ids_with_rna).intersection(set(vcf_samples))
    # subset egan ID and RNA ID links to samples with genotyping data and passing QC 
    wgs_rna_paired_smpls_with_sv_calls_df = wgs_rna_paired_smpls_df[wgs_rna_paired_smpls_df.egan_id.isin(vcf_egan_ids_with_rna)]
    rna_id_pass_qc_sv_calls = wgs_rna_paired_smpls_with_sv_calls_df.rna_id.unique().tolist()
    ge_matrix_flat_chrom_egan_df = pd.merge(ge_mtx_chrom_df, wgs_rna_paired_smpls_with_sv_calls_df, how="inner", on="rna_id")
    # add gene sample pairs 
    ge_matrix_flat_chrom_egan_df["gene_smpl_pair"] = ge_matrix_flat_chrom_egan_df.gene_id + "," + ge_matrix_flat_chrom_egan_df.egan_id

    rna_ids_pass_gt_rna_qc = ge_matrix_flat_chrom_egan_df.rna_id.unique()
    gene_ids_pass_qc = ge_matrix_flat_chrom_egan_df.gene_id.unique()
    print(f"Samples with genotyping data and passing RNA-seq QC: {len(rna_ids_pass_gt_rna_qc)}")
    if ge_matrix_flat_chrom_egan_df.shape[0] != num_misexp_genes * len(rna_ids_pass_gt_rna_qc): 
        raise ValueError("Number of genes and samples passsing QC does not match gene expression matrix size.")

    ### get variants in windows 
    # load variant and gene bed files 
    vrnts_bed = BedTool(vrnts_bed_path)
    genes_bed = BedTool(genes_bed_path)
    # column names for variant window intersection bed file 
    intersect_bed_columns={0:"chrom_gene", 1:"start_gene", 2:"end_gene", 3: "gene_id", 
                           4:"0", 5: "strand", 6:"chrom_vrnt", 7:"start_vrnt", 8:"end_vrnt", 9:"vrnt_id", 10: "AF"}

    vrnts_in_windows_dict = {}
    for direction in ["upstream", "downstream"]: 
        window_size = window_start
        seen_vrnt_gene_pairs = set()
        while window_size <= window_max:
            print(f"{direction}_{window_size}")
            if direction == "upstream": 
                l, r = window_size, 0 
            else: 
                l, r = 0, window_size
            intersect_bed_str = StringIO(str(genes_bed.window(vrnts_bed, l=l, r=r,sw=True)))
            intersect_bed_df = pd.read_csv(intersect_bed_str, sep="\t", header=None).rename(columns=intersect_bed_columns)
            # variant gene ID column 
            intersect_bed_df["vrnt_gene_id"] = intersect_bed_df["vrnt_id"].astype(str) + "," + intersect_bed_df["gene_id"].astype(str)
            # subset to variant gene pairs that have not been seen in previous windows    
            vrnt_gene_pairs_in_window = set(intersect_bed_df.vrnt_gene_id.unique())
            new_vrnt_gene_pairs = vrnt_gene_pairs_in_window - seen_vrnt_gene_pairs
            intersect_bed_new_vrnt_gene_df = intersect_bed_df[intersect_bed_df["vrnt_gene_id"].isin(new_vrnt_gene_pairs)]

            # add variant gene pairs to seen set 
            seen_vrnt_gene_pairs = seen_vrnt_gene_pairs.union(vrnt_gene_pairs_in_window)
            vrnts_in_windows_dict[f"{direction}_{window_size}"] = intersect_bed_new_vrnt_gene_df
            window_size += window_step_size

    # check upstream zero and downstream 0 are the same 
    if not vrnts_in_windows_dict["upstream_0"].equals(vrnts_in_windows_dict["downstream_0"]): 
        print(f"{chrom}: upstream and downstream gene body windows are not identical.")

    intersect_gene_window_str = StringIO(str(genes_bed.window(vrnts_bed, w=gene_window_size)))
    intersect_gene_window_bed_df = pd.read_csv(intersect_gene_window_str, sep="\t", header=None).rename(columns=intersect_bed_columns)
    vrnts_in_windows_dict[f"gene_body_window_{gene_window_size}"] = intersect_gene_window_bed_df

    # create single window for gene body
    vrnts_in_windows_dict["gene_body"] = vrnts_in_windows_dict.pop("upstream_0")
    del vrnts_in_windows_dict["downstream_0"]

    count = 0
    count_carriers = {}
    for z_cutoff in z_cutoff_list: 
        # select misexpression events 
        misexp_df = ge_matrix_flat_chrom_egan_df[(ge_matrix_flat_chrom_egan_df["z-score"] > z_cutoff) & 
                                                          (ge_matrix_flat_chrom_egan_df.TPM > tpm_cutoff)
                                                         ]
        # get gene-sample pairs passing cutoff 
        test_gene_ids = misexp_df.gene_id.unique()
        test_gene_smpl_pairs = misexp_df.gene_smpl_pair.unique()
        # get gene-sample pairs in control group, limit to genes with a misexpression event passing cutoff
        cntrl_gene_smpl_pairs = ge_matrix_flat_chrom_egan_df[~ge_matrix_flat_chrom_egan_df.gene_smpl_pair.isin(test_gene_smpl_pairs) &
                                                             ge_matrix_flat_chrom_egan_df.gene_id.isin(test_gene_ids)
                                                            ].gene_smpl_pair.unique()
        # check overlap between sets is empty 
        if len(set(test_gene_smpl_pairs).intersection(set(cntrl_gene_smpl_pairs))) != 0:
            raise ValueError("Overlap between control and test gene-sample pair sets")
        # check number of test and control gene-pairs matches expected 
        total_test_control_pairs = len(test_gene_smpl_pairs) + len(cntrl_gene_smpl_pairs)
        test_genes_by_smpls = len(test_gene_ids) * len(rna_ids_pass_gt_rna_qc)
        if total_test_control_pairs !=  test_genes_by_smpls: 
            raise ValueError("Number of test and control gene pairs does not match test gene IDs by samples")
        for gene_id in test_gene_ids: 
            print(gene_id)
            test_smpls = set([gene_smpl_pair.split(",")[1] for gene_smpl_pair in test_gene_smpl_pairs if f"{gene_id}," in gene_smpl_pair])
            cntrl_smpls = set([gene_smpl_pair.split(",")[1] for gene_smpl_pair in cntrl_gene_smpl_pairs if f"{gene_id}," in gene_smpl_pair])
            # get genotype info
            gene_id_has_gt_intersect_path = [gene_id_path for gene_id_path in gt_info_root_path.glob(f"{chrom}/*") if gene_id in str(gene_id_path)]
            # check if gene has variants in window
            if len(gene_id_has_gt_intersect_path) == 0:
                for window in vrnts_in_windows_dict.keys():
                    for vrnt_type in ["snp", "indel"]:
                        for maf_range in maf_range_list:  
                            count_carriers[count] = [z_cutoff, gene_id, window, vrnt_type, f"{maf_range[0]}-{maf_range[1]}", 0, len(test_smpls), 0, len(cntrl_smpls)]
                            count += 1 
            elif len(gene_id_has_gt_intersect_path) == 1: 
                gt_info_path = gt_info_root_path.joinpath(f"{chrom}/{chrom}_{gene_id}_vrnts_gts_intersect.tsv")
                gt_info_df = pd.read_csv(gt_info_path, sep="\t")
                for window in vrnts_in_windows_dict.keys():
                    window_df = vrnts_in_windows_dict[window]
                    vrnt_id_in_gene_window = set(window_df[window_df.gene_id == gene_id].vrnt_id.unique())
                    vrnt_id_in_gene_window_df = gt_info_df[gt_info_df.vrnt_id.isin(vrnt_id_in_gene_window)]
                    for vrnt_type in ["snp", "indel"]:
                        vrnt_type_in_gene_window_df = vrnt_id_in_gene_window_df[vrnt_id_in_gene_window_df.vrnt_type == vrnt_type]
                        for maf_range in maf_range_list: 
                            maf_lower, maf_upper = maf_range
                            vrnt_ids_low_af = vrnt_type_in_gene_window_df[(vrnt_type_in_gene_window_df.AF >= maf_lower) &
                                                                        (vrnt_type_in_gene_window_df.AF < maf_upper) 
                                                                       ].vrnt_id.unique()
                            vrnt_ids_high_af = vrnt_type_in_gene_window_df[(vrnt_type_in_gene_window_df.AF <= (1 - maf_lower)) &
                                                                         (vrnt_type_in_gene_window_df.AF > (1 - maf_upper)) 
                                                                        ].vrnt_id.unique()
                            vrnt_ids_in_gene_window = set(vrnt_ids_low_af).union(set(vrnt_ids_high_af))
                            # get list of carriers 
                            vrnt_id_carriers_df = vrnt_id_in_gene_window_df[vrnt_id_in_gene_window_df.vrnt_id.isin(vrnt_ids_in_gene_window)]
                            vrnts_id_carriers = set(vrnt_id_carriers_df.vrnt_id.unique()) 
                            if vrnts_id_carriers != vrnt_ids_in_gene_window: 
                                raise ValueError(f"Missing variant IDs from genotype file: {vrnts_id_carriers - vrnt_ids_in_gene_window}")
                            carriers = vrnt_id_carriers_df.egan_id.unique()
                            # get test and control carriers 
                            test_carriers = test_smpls.intersection(carriers)
                            cntrl_carriers = cntrl_smpls.intersection(carriers)
                            count_carriers[count] = [z_cutoff, gene_id, window, vrnt_type, f"{maf_range[0]}-{maf_range[1]}", len(test_carriers), len(test_smpls), len(cntrl_carriers), len(cntrl_smpls)]
                            count += 1 
            else: 
                raise ValueError(f"{gene_id} has multiple file names: {gene_id_has_gt_intersect_path}")
    
    # write to output
    output_columns=["z_cutoff", "gene_id", "window", "vrnt_type", "maf_range", "misexp_carrier", "misexp_total", "control_carrier", "control_total"]
    count_carriers_df = pd.DataFrame.from_dict(count_carriers, orient="index", columns=output_columns)
    path_out = out_dir_path.joinpath(f"{chrom}_carrier_count.tsv")
    count_carriers_df.to_csv(path_out, sep="\t", index=False)
                
                
if __name__ == "__main__":
    main()
    
    
    