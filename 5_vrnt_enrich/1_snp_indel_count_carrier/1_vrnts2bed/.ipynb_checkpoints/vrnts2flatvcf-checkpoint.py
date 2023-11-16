#!/bin/python
import pandas as pd 
import pysam 
import argparse
from pathlib import Path

def main():
    # arguments passed 
    parser = argparse.ArgumentParser(description="")
    # required 
    parser.add_argument("--chrom", 
                        help ="chromosome", 
                        required = True)
    parser.add_argument("--vcf", 
                        help = "VCF file", 
                        required = True)
    parser.add_argument("--root", 
                        help = "root directory output", 
                        required = True)
    # set arguments to global variable 
    global args
    args = parser.parse_args()
    # load arguments 
    chrom = args.chrom
    vcf_path=args.vcf
    root_dir = args.root 
    
    # check root exists, create if not
    Path(root_dir).mkdir(parents=True, exist_ok=True)
    # variant coordinate bed files 
    flat_vcf_dir = f"{root_dir}/flat_vcf"
    Path(flat_vcf_dir).mkdir(parents=True, exist_ok=True)

    vcf = pysam.VariantFile(vcf_path, mode = "r")
    flat_vcf_path = f"{flat_vcf_dir}/{chrom}_flat_vcf.csv"
    records = vcf.fetch()
    with open(flat_vcf_path, "w") as f: 
        for rec in records: 
            # only keep variants passing MAF cutoff
            af = rec.info["AF"]
            # get chromosome, position, reference and alternative alleles 
            chrom, pos, ref, alt_alleles = rec.chrom, rec.pos, rec.ref, rec.alts
            # check for multiallelic alleles
            if len(alt_alleles) > 1: 
                raise ValueError(f"Multiallelic entry in vcf at: {chrom}, {pos}")
            else: 
                alt = alt_alleles[0]
            vcf_vrnt_id = f"{chrom}:{pos}:{ref}:{alt}"
            # check ref and alt alleles do not contain N or . 
            if any(b not in "ATGC" for b in ref+alt): 
                raise ValueError(f"Variant {vcf_variant_id} contains unassigned nucleotides.")
            gts = [s["GT"] for s in rec.samples.values()]
            [f.write(f"{','.join(str(x) for x in [vcf_vrnt_id, vcf.header.samples[i], gt, af])}\n") for i, gt in enumerate(gts)]
        

if __name__ == "__main__":
    main()