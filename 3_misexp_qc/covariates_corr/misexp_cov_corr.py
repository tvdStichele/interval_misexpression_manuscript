#!/bin/python
import pandas as pd 
from pathlib import Path
import scipy.stats as stats
import argparse

def main():
    # arguments passed 
    parser = argparse.ArgumentParser(description="")
    # required 
    parser.add_argument("--name", 
                        help ="Name of dataset", 
                        required = True)
    parser.add_argument("--expression", 
                        help ="Expression matrix of inactive genes", 
                        required = True)
    parser.add_argument("--covariates", 
                        help ="GTEx eqtl covariates", 
                        required = True)
    parser.add_argument("--output", 
                        help ="Output directory", 
                        required = True)
    # set arguments to global variable 
    global args
    args = parser.parse_args()
    
    # load arguments 
    name = args.name
    tpm_inactive_path = args.expression
    covariates_path = args.covariates
    out_dir = args.output
    
    # load gene expression 
    tpm_inactive_df = pd.read_csv(tpm_inactive_path, sep="\t")
    inactive_gene_ids = tpm_inactive_df.gene_id.unique()
    tpm_inactive_df = tpm_inactive_df.set_index("gene_id")
    # reformat TPM matrix 
    tpm_inactive_tp_df = tpm_inactive_df.transpose()
    tpm_inactive_tp_df = tpm_inactive_tp_df.reset_index(drop=False).rename(columns={"index":"rna_id"})

    # add covariates 
    covariates_df = pd.read_csv(covariates_path, sep="\t")
    covariates = [col for col in covariates_df.columns if col not in ["smpl_id", "rna_id"]]

    tpm_inactive_covariates_df = pd.merge(tpm_inactive_tp_df, 
                                          covariates_df, 
                                          on="rna_id", 
                                          how="inner")

    # calculate spearman correlation between covariates and gene expression 
    gene_cov_corr_dict, count = {}, 0
    for gene_id in inactive_gene_ids: 
        for covariate in covariates: 
            correlation, pval = stats.spearmanr(tpm_inactive_covariates_df[gene_id], 
                                                tpm_inactive_covariates_df[covariate], 
                                                nan_policy="omit"
                                               )
            gene_cov_corr_dict[count] = [gene_id, covariate, correlation, pval, name]
            count += 1
    # generate dataframe 
    columns=["gene_id", "covariate", "spearman", "pval", "name"]
    gene_cov_corr_df = pd.DataFrame.from_dict(gene_cov_corr_dict, orient="index", columns=columns)
    # output 
    out_path = Path(out_dir)
    gene_cov_corr_path = out_path.joinpath(f"{name}_gene_cov_corr.tsv")
    gene_cov_corr_df.to_csv(gene_cov_corr_path, sep="\t", index=False)
    
    
if __name__ == "__main__":
    main()