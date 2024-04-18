#!/bin/python
import pandas as pd 
import numpy as np
from pathlib import Path
import argparse

def main():
    # arguments passed 
    parser = argparse.ArgumentParser(description="")
    # required 
    parser.add_argument("--chrom", 
                        help ="chromosome", 
                        required = True)
    parser.add_argument("--sv_info", 
                        help = "structural variants information file", 
                        required = True)
    parser.add_argument("--vep_all", 
                        help = "File with variant all VEP consequences", 
                        required = True)
    parser.add_argument("--vep_msc", 
                        help = "File with VEP MSC consequences", 
                        required = True)
    parser.add_argument('--z_cutoff', 
                        nargs='+', 
                        type=float, 
                        required = True,
                        help = "list of TPM cutoffs")
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
    chrom = args.chrom
    sv_info_path = args.sv_info
    vep_all_path = args.vep_all
    vep_msc_path = args.vep_msc
    z_cutoff_list = args.z_cutoff
    af_cutoff_list = args.af_range
    root_dir = args.root
    af_cutoff_bins = [(af[1], af_cutoff_list[af[0]+1]) for af in enumerate(af_cutoff_list[:-1])]
                       
    print("Inputs:")
    print(f"- Chromosome: {chrom}")
    print(f"- AF bins: {af_cutoff_bins}")
    print(f"- Z-score cutoffs: {z_cutoff_list}")
    print("")
    
    vep_msc_ranks = [
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "transcript_amplification",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        "splice_donor_5th_base_variant", 
        "splice_region_variant",
        "splice_donor_region_variant", 
        "splice_polypyrimidine_tract_variant"
        "incomplete_terminal_codon_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "synonymous_variant",
        "coding_sequence_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "intron_variant",
        "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TFBS_ablation",
        "TFBS_amplification",
        "TF_binding_site_variant",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "feature_elongation",
        "regulatory_region_variant",
        "feature_truncation",
        "intergenic_variant",
        "sequence_variant",
        "no_predicted_effect"]
    
    # from all VEP consequences that have no gene annotation
    vep_msc_not_linked_to_gene = ['TFBS_ablation', 
                                  'TF_binding_site_variant',
                                  'regulatory_region_variant', 
                                  'TFBS_amplification',
                                  'intergenic_variant', 
                                  'regulatory_region_ablation',
                                  'regulatory_region_amplification']

    root_dir_path = Path(root_dir)
    tpm_cutoff = 0.5

    # load gene expression matrix 
    ge_matrix_flat_chrom_egan_path = root_dir_path.joinpath(f"express_mtx/{chrom}_ge_mtx_flat.tsv")
    ge_matrix_flat_chrom_egan_df = pd.read_csv(ge_matrix_flat_chrom_egan_path, sep="\t")
    # add gene-sample pairs
    ge_matrix_flat_chrom_egan_df["gene_smpl_pair"] = ge_matrix_flat_chrom_egan_df.gene_id + "," + ge_matrix_flat_chrom_egan_df.rna_id
    # carrier information 
    sv_intersect_express_info_path = root_dir_path.joinpath(f"express_carrier_info/{chrom}_express_carrier_info.tsv")
    sv_intersect_express_info_df = pd.read_csv(sv_intersect_express_info_path, sep="\t")
    # add gene-sample pairs  
    sv_intersect_express_info_df["gene_smpls_id"] = sv_intersect_express_info_df.gene_id.astype(str) + "," + sv_intersect_express_info_df.rna_id.astype(str)
    # RNA IDs with EGAN ID match 
    sv_intersect_express_info_path = root_dir_path.joinpath("egan_rna_smpls/egan_rna_ids_paired_pass_qc.tsv")
    rna_id_pass_qc_sv_calls = pd.read_csv(sv_intersect_express_info_path, sep="\t").rna_id.unique()
    # SV types 
    sv_info_df = pd.read_csv(sv_info_path, sep="\t", dtype={"plinkID":str})
    sv_info_id_af_df = sv_info_df[["plinkID", "AF", "SVTYPE"]].rename(columns={"plinkID":"vrnt_id"})
    sv_types_list = sv_info_df.SVTYPE.unique().tolist()
    
    # annotate variant-gene pairs with VEP consequence 
    vep_all_consq_df = pd.read_csv(vep_all_path, sep="\t")
    # variant gene pairs 
    vrnt_gene_pairs_df = sv_intersect_express_info_df[["vrnt_id", "gene_id"]].drop_duplicates()
    # collapse to gene-variant consequence pairs 
    vep_all_consq_drop_dups_df = vep_all_consq_df[["vrnt_id", "gene_id", "Consequence"]].drop_duplicates()
    # add variant-gene consequences 
    vrnt_gene_pairs_consq_df = pd.merge(vrnt_gene_pairs_df, 
                                        vep_all_consq_drop_dups_df, 
                                        on=["vrnt_id", "gene_id"], 
                                        how="left").fillna("no_predicted_effect")
    # group consequences together for each variant-gene pair
    vrnt_gene_pairs_consq_df["gene_consequence"] = vrnt_gene_pairs_consq_df.groupby(["vrnt_id", "gene_id"])['Consequence'].transform(lambda x: ','.join(x))
    # remove duplicates arising from variant have different consequences across gene transcripts 
    vrnt_gene_pair_consq_collapse_df = vrnt_gene_pairs_consq_df.drop(columns=["Consequence"]).drop_duplicates()
    if vrnt_gene_pair_consq_collapse_df.shape[0] != vrnt_gene_pairs_df.shape[0]: 
        raise ValueError("Variant gene pair number does not match number of variant gene pairs with consequence.")
    # assign each variant-gene pair a unique consequence based on VEP rank  
    vrnt_gene_pair_msc_consq = []
    for index, row in vrnt_gene_pair_consq_collapse_df.iterrows(): 
        gene_consequence = row["gene_consequence"]
        gene_consequence_list = gene_consequence.split(",")
        for consq in vep_msc_ranks: 
            if consq in gene_consequence_list: 
                break 
        vrnt_gene_pair_msc_consq.append(consq)
    if vrnt_gene_pair_consq_collapse_df.shape[0] != len(vrnt_gene_pair_msc_consq): 
        raise ValueError("Number of MSC per gene does not match number of variant-gene pairs.")
    vrnt_gene_pair_consq_collapse_df["consequence"] = vrnt_gene_pair_msc_consq
    vrnt_gene_pair_consq_collapse_df = vrnt_gene_pair_consq_collapse_df.drop(columns=["gene_consequence"])
    # split into variants with annotated gene effect and no predicted gene effect
    vrnt_gene_effect_df = vrnt_gene_pair_consq_collapse_df[vrnt_gene_pair_consq_collapse_df.consequence != "no_predicted_effect"]
    no_predicted_effect_df = vrnt_gene_pair_consq_collapse_df[vrnt_gene_pair_consq_collapse_df.consequence == "no_predicted_effect"]
    # load VEP MSC annotations
    vep_msc_df = pd.read_csv(vep_msc_path, sep="\t").rename(columns={"Uploaded_variation": "vrnt_id", "Consequence": "msc"})
    # merge VEP MSC with no predicted effect 
    no_prediced_effect_df = pd.merge(no_predicted_effect_df[["vrnt_id", "gene_id"]], 
                                     vep_msc_df[["vrnt_id", "msc"]], 
                                     on="vrnt_id", 
                                     how="left")
    # check for NaNs 
    if no_prediced_effect_df.msc.isnull().values.any(): 
        raise ValueError("Variants missing VEP most severe consequence.")
    no_prediced_effect_df["consequence"] = np.where(no_prediced_effect_df.msc.isin(vep_msc_not_linked_to_gene), 
                                                    no_prediced_effect_df.msc, 
                                                    "no_predicted_effect")
    no_prediced_effect_trunc_df = no_prediced_effect_df[["vrnt_id", "gene_id", "consequence"]]
    # combine variants with annotated gene effect and variants with no predicted effect with updated regulatory consequence
    vrnt_gene_pair_consq_msc_added_df = pd.concat([vrnt_gene_effect_df, no_prediced_effect_trunc_df])
    if vrnt_gene_pair_consq_msc_added_df[["vrnt_id", "gene_id"]].drop_duplicates().shape[0] != vrnt_gene_pairs_df.shape[0]: 
        raise ValueError("Number of gene pairs in expression input does not match number of gene pairs with annotated consequence.")
    # add consequence to carrier information 
    sv_intersect_express_info_gene_msc_df = pd.merge(sv_intersect_express_info_df, 
                                                     vrnt_gene_pair_consq_msc_added_df, 
                                                     on=["vrnt_id", "gene_id"], 
                                                     how="left")
    # check number of rows not changed 
    if sv_intersect_express_info_gene_msc_df.shape[0] != sv_intersect_express_info_df.shape[0]: 
        raise ValueError("Added entries to sample-gene-variant expression dataframe.")
    # check for nans 
    if sv_intersect_express_info_gene_msc_df.consequence.isnull().values.any(): 
        raise ValueError("Missing consequence annotations for some variant-gene pairs.")
    ### count carriers in test and controls 
    carrier_count = {}
    gene_smpl_pair_carrier = {}
    gene_smpl_count = 0 
    for af_range in af_cutoff_bins:
        af_lower, af_upper = af_range
        for z_cutoff in z_cutoff_list: 
            carrier_count[gene_smpl_count] = [chrom, z_cutoff, f"{round(af_lower*100)}-{round(af_upper*100)}"]
            ### get control and test gene-sample pairs based on misexpression cutoff 
            # get gene-sample pairs passing z-score cutoff and misexpression cutoff
            z_cutoff_df = ge_matrix_flat_chrom_egan_df[(ge_matrix_flat_chrom_egan_df["z-score"] > z_cutoff) & 
                                                       (ge_matrix_flat_chrom_egan_df["TPM"] > tpm_cutoff)
                                                      ]
            test_gene_ids = z_cutoff_df.gene_id.unique()
            test_gene_smpl_pairs = z_cutoff_df.gene_smpl_pair.unique()
            # get gene-sample pairs in control group 
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
            carrier_count[gene_smpl_count] += [len(test_gene_ids), len(rna_id_pass_qc_sv_calls), len(test_gene_smpl_pairs), len(cntrl_gene_smpl_pairs)]

            ### count carriers in control and test groups 
            # get misexpression carriers 
            misexpress_carriers_df = sv_intersect_express_info_gene_msc_df[(sv_intersect_express_info_gene_msc_df.gene_smpls_id.isin(test_gene_smpl_pairs)) & 
                                                                  (sv_intersect_express_info_gene_msc_df.AF >= af_lower) & 
                                                                  (sv_intersect_express_info_gene_msc_df.AF < af_upper) &
                                                                  (sv_intersect_express_info_gene_msc_df.genotype.isin(['(0, 1)', '(1, 1)']))
                                                                 ].copy()
            # get non-misexpression carriers
            cntrl_carriers_df = sv_intersect_express_info_gene_msc_df[(sv_intersect_express_info_gene_msc_df.gene_smpls_id.isin(cntrl_gene_smpl_pairs)) &
                                                                      (sv_intersect_express_info_gene_msc_df.AF >= af_lower) & 
                                                                      (sv_intersect_express_info_gene_msc_df.AF < af_upper) &
                                                                      (sv_intersect_express_info_gene_msc_df.genotype.isin(['(0, 1)', '(1, 1)']))
                                                                     ].copy()
            # add number of carriers 
            gene_smpl_misexpress_carriers = misexpress_carriers_df.gene_smpls_id.unique()
            gene_smpl_cntrl_carriers = cntrl_carriers_df.gene_smpls_id.unique()
            carrier_count[gene_smpl_count] += [len(gene_smpl_misexpress_carriers), len(gene_smpl_cntrl_carriers)]
            # check overlap between sets is empty 
            if len(set(gene_smpl_misexpress_carriers).intersection(set(gene_smpl_cntrl_carriers))) != 0:
                raise ValueError("Overlap between control and test gene-sample pair carrier sets")

            ### count carriers for SV types
            for sv_type in sv_types_list:
                misexpress_carriers_sv_type_df = misexpress_carriers_df[misexpress_carriers_df.SVTYPE == sv_type]
                cntrl_carriers_sv_type_df = cntrl_carriers_df[cntrl_carriers_df.SVTYPE == sv_type]
                sv_type_misexp_carriers = len(misexpress_carriers_sv_type_df.gene_smpls_id.unique())
                sv_type_cntrl_carriers = len(cntrl_carriers_sv_type_df.gene_smpls_id.unique())
                carrier_count[gene_smpl_count] += [sv_type_misexp_carriers, sv_type_cntrl_carriers]
                for msc in vep_msc_ranks: 
                    msc_type_misexp_carriers = len(misexpress_carriers_sv_type_df[misexpress_carriers_sv_type_df.consequence == msc].gene_smpls_id.unique())
                    msc_type_cntrl_carriers = len(cntrl_carriers_sv_type_df[cntrl_carriers_sv_type_df.consequence == msc].gene_smpls_id.unique())
                    carrier_count[gene_smpl_count] += [msc_type_misexp_carriers, msc_type_cntrl_carriers]
            gene_smpl_count += 1 

    ### write results to file 
    vrnt_type_cols = []
    for sv_type in sv_types_list: 
        vrnt_type_cols += [f"{sv_type}_misexp", f"{sv_type}_contrl"]
        for msc in vep_msc_ranks:
            vrnt_type_cols += [f"{sv_type}_{msc}_misexp", f"{sv_type}_{msc}_contrl"]

    columns = ["chrom", "z_cutoff", "maf_bin", "misexp_genes", "smpls_pass_qc", "total_misexp", "total_control", "all_sv_misexp", "all_sv_contrl"] + vrnt_type_cols
    carrier_count_df = pd.DataFrame.from_dict(carrier_count, orient="index", columns=columns)
    # write gene sample pair carrier counts 
    carrier_count_dir = root_dir_path.joinpath("carrier_count_gene_msc_reg_af50")
    Path(carrier_count_dir).mkdir(parents=True, exist_ok=True)
    carrier_count_path = carrier_count_dir.joinpath(f"{chrom}_carrier_count_gene_msc.tsv")
    carrier_count_df.to_csv(carrier_count_path, sep="\t", index=False)
    
    
if __name__ == "__main__":
    main()