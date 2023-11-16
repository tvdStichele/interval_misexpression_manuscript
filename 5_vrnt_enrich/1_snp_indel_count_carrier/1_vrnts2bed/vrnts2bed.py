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
                        help = "VCF file. Multiallelic entries on separate lines.", 
                        required = True)
    parser.add_argument("--max_indel_length", 
                        help = "Maximum indel length", 
                        type = int,
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
    max_indel_length=args.max_indel_length
    root_dir = args.root 
    
    # check root exists, create if not, create variant coordinate bed files 
    root_dir_path = Path(root_dir)
    vrnts_bed_dir = root_dir_path.joinpath("vrnts_bed")
    vrnts_bed_dir.mkdir(parents=True, exist_ok=True)

    vcf = pysam.VariantFile(vcf_path, mode = "r")
    vrnts_bed_path = vrnts_bed_dir.joinpath(f"{chrom}_vrnts.bed")
    records = vcf.fetch()
    with open(vrnts_bed_path, "w") as f: 
        for rec in records: 
            af = rec.info["AF"]
            # get chromosome, position, reference and alternative alleles 
            chrom, pos, ref, alt_alleles = rec.chrom, rec.pos, rec.ref, rec.alts
            chrom_num = chrom.split("chr")[1]
            # check for multiallelic alleles
            if len(alt_alleles) > 1: 
                raise ValueError(f"Multiallelic entry in vcf at: {chrom}, {pos}")
            else: 
                alt = alt_alleles[0]
            vcf_vrnt_id = f"{rec.chrom}:{rec.pos}:{ref}:{alt}"
            # check ref and alt alleles do not contain N or . 
            if any(b not in "ATGC" for b in ref+alt): 
                raise ValueError(f"Variant {vcf_vrnt_id} contains unassigned nucleotides.")
            # convert to bed 
            # deletions (vcf2bed includes base preceding deletion)
            if len(ref) > len(alt):
                if len(alt) > 1: 
                    raise ValueError(f"Variant {vcf_vrnt_id} alt allele length is greater than one: {len(alt)}") 
                start = pos - 1 
                end = pos+len(ref) - 1
                length = len(ref) - len(alt)
            # snvs and insertions (vcf2bed yields a one-base bed element)
            else: 
                start = pos - 1 
                end = pos
                length = len(alt) - len(ref)
            # check length 
            if length < 0: 
                raise ValueError(f"Negative variant length: {vcf_vrnt_id}")
            # remove indels longer than max length  
            elif length > max_indel_length: 
                continue 
            else: 
                line = "\t".join(str(x) for x in [chrom_num, start, end, vcf_vrnt_id, af])
                f.write(f"{line}\n")


if __name__ == "__main__":
    main()