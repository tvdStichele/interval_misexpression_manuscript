#!/bin/python
import pandas as pd 
import pysam 
import argparse
from pathlib import Path


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
    parser.add_argument("--gencode", 
                        help = "gencode file", 
                        required = True)
    parser.add_argument("--root", 
                        help = "root directory output", 
                        required = True)
    # constants 
    tpm_cutoff = 0.5
    z_score_cutoff = 2
    # set arguments to global variable 
    global args
    args = parser.parse_args()
    # load arguments 
    ge_matrix_flat_path = args.express
    chrom = args.chrom
    gencode_path=args.gencode
    root_dir = args.root 
    
    ### write genes to bed file 
    # check root exists - if not make 
    root_dir_path =  Path(root_dir)
    root_dir_path.mkdir(parents=True, exist_ok=True)
    # gene coordinate bed files 
    gene_bed_dir = root_dir_path.joinpath("genes_bed")
    Path(gene_bed_dir).mkdir(parents=True, exist_ok=True)
    # load gene expression data
    ge_matrix_flat_df = pd.read_csv(ge_matrix_flat_path, sep=",")
    misexp_df = ge_matrix_flat_df[(ge_matrix_flat_df.TPM > tpm_cutoff) & 
                                  (ge_matrix_flat_df["z-score"] > z_score_cutoff)]                          
    gene_id_pass_qc_set = set(misexp_df.gene_id.unique())
    print(f"Gene IDs passing QC with at least one misexpression event: {len(gene_id_pass_qc_set)}")
    
    # write output to bed
    gene_bed_path = gene_bed_dir.joinpath(f"{chrom}_genes.bed")
    gene_id_in_bed = []
    # load gencode gtf
    with open(gene_bed_path, "w") as bed_out:
        for gtf in pysam.TabixFile(gencode_path).fetch(chrom, parser = pysam.asGTF()):
            if gtf.feature == "gene" and gtf.gene_id.split('.')[0] in gene_id_pass_qc_set:
                # collect all entries with gene ID 
                gene_id_list, chrom_list, start_list, end_list, strand_list = [], [], [], [], []
                gene_id_list.append(gtf.gene_id.split('.')[0])
                chrom_list.append(gtf.contig)
                start_list.append(gtf.start)
                end_list.append((gtf.end))
                strand_list.append(gtf.strand)
                # check or write output
                if len(gene_id_list) > 1 or len(chrom_list) > 1 or len(start_list) > 1 or len(end_list) > 1:
                    print(f"{gene_id} has multiple entries in gencode file - not included in output file.")
                elif len(chrom_list) == 0 or len(start_list) == 0 or len(end_list) == 0: 
                    print(f"{gene_id} has no entries in gencode file - not included in output file.")
                else: 
                    gene_id, gtf_chrom, start, end, strand = gene_id_list[0], chrom_list[0], int(start_list[0]) - 1, int(end_list[0]), strand_list[0]
                    gene_id_in_bed.append(gene_id)
                    chrom_num = gtf_chrom.split("chr")[1]
                    # write gene to bed file 
                    bed_out.write(f"{chrom_num}\t{start-1}\t{end}\t{gene_id}\t0\t{strand}\n")
    bed_out.close()
    print(f"Number of genes pass QC on {chrom}: {len(gene_id_in_bed)}")
    
    
if __name__ == "__main__":
    main()