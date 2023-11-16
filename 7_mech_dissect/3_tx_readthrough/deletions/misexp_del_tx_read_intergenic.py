#!/bin/python
import pandas as pd
from pathlib import Path
import pysam
from pybedtools import BedTool
from io import StringIO
import argparse

def main():
    # arguments passed 
    parser = argparse.ArgumentParser(description="")
    # required
    parser.add_argument("--vrnt_id", 
                        help ="Variant ID", 
                        required = True)
    parser.add_argument("--tx_read_info", 
                        help ="Information on misexpression-associated variants", 
                        required = True)
    parser.add_argument("--bam_dir", 
                        help ="Bams directory", 
                        required = True)
    parser.add_argument("--out", 
                    help ="Output", 
                    required = True)
    
    ### inputs 
    # set arguments to global variable 
    global args
    args = parser.parse_args()
    vrnt_id=args.vrnt_id
    tx_read_info_path = args.tx_read_info
    misexp_cntrl_bams_dir = args.bam_dir
    out_dir = args.out
    
    misexp_cntrl_bams_path = Path(misexp_cntrl_bams_dir)
    misexp_tx_read_vrnt_feat_df = pd.read_csv(tx_read_info_path, sep="\t", dtype={"vrnt_id": str})

    intergenic_cov_df_list = []

    # get SV information
    misexp_vrnt_id_df = misexp_tx_read_vrnt_feat_df[misexp_tx_read_vrnt_feat_df.vrnt_id == vrnt_id]
    sv_start = misexp_vrnt_id_df.sv_start.unique()[0]
    sv_end = misexp_vrnt_id_df.sv_end.unique()[0]
    chrom = misexp_vrnt_id_df.chrom.unique()[0]
    seq_name = chrom.split("chr")[1]
    # misexpressed genes 
    misexp_genes = misexp_vrnt_id_df[misexp_vrnt_id_df.vrnt_id == vrnt_id].gene_id.unique()
    for gene_id in misexp_genes:
        # get gene information
        misexp_vrnt_id_gene_df = misexp_vrnt_id_df[misexp_vrnt_id_df.gene_id == gene_id]
        gene_strand = misexp_vrnt_id_gene_df.gene_strand.unique()[0]
        gene_start = misexp_vrnt_id_gene_df.gene_start.unique()[0]
        gene_end = misexp_vrnt_id_gene_df.gene_end.unique()[0]
        misexp_rna_ids = misexp_vrnt_id_gene_df.rna_id.unique()
        # intergenic region
        if gene_strand == '+': 
            start = sv_end
            end = gene_start
        elif gene_strand == "-":
            start = gene_end
            end = sv_start
        else:
            raise ValueError(f"{gene_strand} not recognised.") 
        # intergenic region bed
        vrnt_gene_id = f"{vrnt_id}_{gene_id}"
        intergenic_reg = BedTool(f'{chrom} {start} {end} {vrnt_gene_id} 0 {gene_strand}', from_string=True)
        # variant bams 
        vrnt_bams_dir = misexp_cntrl_bams_path.joinpath(f"{vrnt_id}")
        vrnt_bams_path_list = vrnt_bams_dir.glob("*.bam")
        # compute coverage metrics over intergenic regions
        for path in vrnt_bams_path_list:  
            bam_bed = BedTool(path)
            # compute coverage metrics only on strand with misexpressed gene 
            intergenic_cov_str = StringIO(str(intergenic_reg.coverage(bam_bed, split=True, s=True)))
            intergenic_cov_df = pd.read_csv(intergenic_cov_str, sep="\t", header=None)
            if intergenic_cov_df.shape[0] != 1: 
                raise ValueError("Empty intergenic region")
            intergenic_cov_df["vrnt_id"] = vrnt_id
            intergenic_cov_df["gene_id"] = gene_id
            rna_id = path.name.split(".")[0]
            intergenic_cov_df["rna_id"] = rna_id
            intergenic_cov_df_list.append(intergenic_cov_df)

    intergenic_reads_columns = {0: "chrom", 1:"sv_start", 2:"sv_end", 3:"vrnt_gene_id", 
                                4:"score", 5:"misexp_gene_strand",6:"features", 7:"coverage", 
                                8:"total_len", 9: "cov_fraction"}
    tx_read_vrnts_intergenic_cov_df = pd.concat(intergenic_cov_df_list)
    tx_read_vrnts_intergenic_cov_df = tx_read_vrnts_intergenic_cov_df.rename(columns=intergenic_reads_columns)
    # output directory
    out_dir_path = Path(out_dir)
    vrnt_out_path = out_dir_path.joinpath(f"{vrnt_id}")
    vrnt_out_path.mkdir(parents=True, exist_ok=True)
    # write to file 
    tx_read_vrnts_intergenic_cov_path = vrnt_out_path.joinpath(f"{vrnt_id}_intergenic_cov.tsv")
    tx_read_vrnts_intergenic_cov_df.to_csv(tx_read_vrnts_intergenic_cov_path, sep="\t", index=False)
    
    
if __name__ == "__main__":
    main()