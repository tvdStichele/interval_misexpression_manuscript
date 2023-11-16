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
    parser.add_argument("--chrom", 
                        help ="chromosome", 
                        required = True)
    parser.add_argument("--gene_bed", 
                        help = "Bed file with genes to test", 
                        required = True)
    parser.add_argument("--vrnts_bed", 
                        help = "Bed file with variants in VCF", 
                        required = True)
    parser.add_argument("--window", 
                        help = "Size of window around gene", 
                        type=int,
                        required = True)
    parser.add_argument("--root", 
                        help = "root directory output", 
                        required = True)
    # set arguments to global variable 
    global args
    args = parser.parse_args()
    # load arguments 
    chrom = args.chrom
    vrnts_bed_path = args.vrnts_bed
    genes_bed_path = args.gene_bed
    window = args.window
    root_dir = args.root 
    
    # constants 
    window_intersect_cols={0:"chrom_gene", 1:"start_gene", 2:"end_gene", 3: "gene_id", 
                           4:"score", 5:"strand", 6:"chrom_vrnt", 7:"start_vrnt", 8:"end_vrnt", 9:"vrnt_id", 10:"AF"}

    # check root exists
    root_dir_path = Path(root_dir)
    root_dir_path.mkdir(parents=True, exist_ok=True)

    ### Find variants intersecting gene windows 
    intersect_bed_dir = root_dir_path.joinpath("intersect_beds")
    Path(intersect_bed_dir).mkdir(parents=True, exist_ok=True)
    intersect_bed_path = intersect_bed_dir.joinpath(f"{chrom}_gene_vrnts_intersect.tsv")
    vrnts_bed = BedTool(vrnts_bed_path)
    gene_bed= BedTool(genes_bed_path)
    print(f"Finding variant IDs in gene windows ...")
    window_intersect_vrnts_str = StringIO(str(gene_bed.window(vrnts_bed, w=window)))
    window_intersect_vrnts_df = pd.read_csv(window_intersect_vrnts_str, sep="\t", header=None).rename(columns=window_intersect_cols)
    window_intersect_vrnts_df.to_csv(intersect_bed_path, sep="\t", index=False)
    print(f"Completed.")
    # genes with variants in windows 
    intersect_genes_dir = root_dir_path.joinpath("intersect_genes")
    Path(intersect_genes_dir).mkdir(parents=True, exist_ok=True)
    intersect_genes_path = intersect_genes_dir.joinpath(f"{chrom}_gene_ids_intersect.txt")
    intersect_gene_ids_list = window_intersect_vrnts_df.gene_id.unique()
    print(f"{chrom} number of genes with variant in windows: {len(intersect_gene_ids_list)}")
    # write gene list 
    print("Writing genes to file ...")
    with open(intersect_genes_path, 'w') as genes_out:
        for gene_id in intersect_gene_ids_list: 
            genes_out.write(f"{gene_id}\n")
    print("Completed.")
    
    
if __name__ == "__main__":
    main()
