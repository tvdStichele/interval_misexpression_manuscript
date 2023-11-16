#!/bin/python
import pandas as pd 
from pathlib import Path
import pysam
import argparse

SWAP_SAMPLES_DICT = {"INT_RNA7879032": "INT_RNA7879033", 
                    "INT_RNA7879033": "INT_RNA7879032", 
                    "INT_RNA7960192": "INT_RNA7960193", 
                    "INT_RNA7960193": "INT_RNA7960192",
                    "INT_RNA7709692": "INT_RNA7709693",
                    "INT_RNA7709693": "INT_RNA7709692",
                    "INT_RNA7710161": "INT_RNA7710162", 
                    "INT_RNA7710162": "INT_RNA7710161",
                    "INT_RNA7710163": "INT_RNA7710164", 
                    "INT_RNA7710164": "INT_RNA7710163"}


CHROMOSOMES_NUM = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
                   '13', '14', '15', '16', '17', '18','19', '20', '21', '22', 'X', 'Y']

def main():
    # arguments passed 
    parser = argparse.ArgumentParser(description="")
    # required
    parser.add_argument("--vrnt_id", 
                        help ="Variant ID", 
                        required = True)
    parser.add_argument("--star_dir", 
                    help ="STAR 2nd pass output directory", 
                    required = True)
    parser.add_argument("--rna_pass", 
                    help ="RNA IDs passing QC", 
                    required = True)
    parser.add_argument("--misexp_feat", 
                        help ="Information on misexpression-associated variants", 
                        required = True)
    parser.add_argument("--out", 
                    help ="Output", 
                    required = True)
    
    ### inputs 
    # set arguments to global variable 
    global args
    args = parser.parse_args()
    vrnt_id=args.vrnt_id
    star_pass_2_out = args.star_dir
    rna_ids_pass_qc_path = args.rna_pass
    misexp_feat_path = args.misexp_feat
    out_dir = args.out
    
    window = 1000000 # cis window set to 1 Mb

    # RNA IDs passing QC 
    rna_ids_pass_qc = pd.read_csv(rna_ids_pass_qc_path, sep="\t", header=None)[0].tolist()
    
    # output directory
    out_dir_path = Path(out_dir)
    # vrnt info 
    misexp_vrnt_feat_df = pd.read_csv(misexp_feat_path, sep="\t")
    # STAR output 
    star_pass_2_path = Path(star_pass_2_out)
    # get SV coordinates
    misexp_vrnt_id_df = misexp_vrnt_feat_df[["vrnt_id", "chrom", "sv_start", "sv_end"]].drop_duplicates()
    vrnt_id_df = misexp_vrnt_id_df[misexp_vrnt_id_df.vrnt_id == vrnt_id]
    # iterate over variants 
    for index, row in vrnt_id_df.iterrows(): 
        vrnt_id = row["vrnt_id"]
        chrom = row["chrom"]
        seq_name = chrom.split("chr")[1]
        sv_start = row["sv_start"]
        sv_end = row["sv_end"]
        # cis window 
        start = sv_start - window 
        start = start if start > 0 else 0
        end = sv_end + window
        # output bam directory 
        bam_out_dir_path = out_dir_path.joinpath(vrnt_id)
        bam_out_dir_path.mkdir(parents=True, exist_ok=True)
        # get misexpressed samples 
        misexp_rna_ids = misexp_vrnt_feat_df[misexp_vrnt_feat_df.vrnt_id == vrnt_id].rna_id.unique()
        for rna_id in rna_ids_pass_qc: 
            if rna_id in misexp_rna_ids: 
                smpl_type = "misexp"
            else: 
                smpl_type = "cntrl"

            # swap RNA ID sample if in swaps 
            if rna_id in SWAP_SAMPLES_DICT.keys(): 
                bam_rna_id = SWAP_SAMPLES_DICT[rna_id]
            else: 
                bam_rna_id = rna_id

            # subset bam file
            bam_path = star_pass_2_path.joinpath(f"{bam_rna_id}/{bam_rna_id}.Aligned.sortedByCoord.out.bam")

            # load input .bam
            bamfile = pysam.AlignmentFile(bam_path, "rb")
            # load output .bam - correct the sample swap
            output_bamfile_path = bam_out_dir_path.joinpath(f"{rna_id}.{smpl_type}.{vrnt_id}.chr{seq_name}.{start}.{end}.Aligned.bam")
            output = pysam.AlignmentFile(output_bamfile_path, "wb", template=bamfile)

            # subset bam file to specific location 
            bamfile_reads = bamfile.fetch(seq_name, start, end)
            for read in bamfile_reads:
                output.write(read)
            bamfile.close()
            output.close()

            # index and sort bam 
            output_bamfile_srtd_path = bam_out_dir_path.joinpath(f"{rna_id}.{smpl_type}.{vrnt_id}.chr{seq_name}.{start}.{end}.Aligned.sorted.bam")
            pysam.sort("-o", str(output_bamfile_srtd_path), str(output_bamfile_path))
            pysam.index(str(output_bamfile_srtd_path))

            # adjust chromosome name in header 
            full_output_path = bam_out_dir_path.joinpath(f"{rna_id}.{smpl_type}.{vrnt_id}.chr{seq_name}.{start}.{end}.Aligned.sorted.chr_names.bam")
            prefix = 'chr'
            output_bamfile =  pysam.AlignmentFile(output_bamfile_srtd_path, "rb")
            new_head = output_bamfile.header.to_dict()
            for seq in new_head['SQ']:
                if seq["SN"] in CHROMOSOMES_NUM:
                    seq['SN'] = prefix + seq['SN']

            # create output BAM with newly defined header
            with pysam.AlignmentFile(full_output_path, "wb", header=new_head) as outf:
                for read in output_bamfile.fetch():
                    prefixed_chrom = prefix + read.reference_name
                    a = pysam.AlignedSegment(outf.header)
                    a.query_name = read.query_name
                    a.query_sequence = read.query_sequence
                    a.reference_name = prefixed_chrom
                    a.flag = read.flag
                    a.reference_start = read.reference_start
                    a.mapping_quality = read.mapping_quality
                    a.cigar = read.cigar
                    a.next_reference_id = read.next_reference_id
                    a.next_reference_start = read.next_reference_start
                    a.template_length = read.template_length
                    a.query_qualities = read.query_qualities
                    a.tags = read.tags
                    outf.write(a)
            # index output with chromosome names 
            pysam.index(str(full_output_path))

            # remove intermediate files 
            Path.unlink(output_bamfile_path)
            Path.unlink(output_bamfile_srtd_path)
            output_bamfile_index_path = bam_out_dir_path.joinpath(f"{rna_id}.{smpl_type}.{vrnt_id}.chr{seq_name}.{start}.{end}.Aligned.sorted.bam.bai")
            Path.unlink(output_bamfile_index_path)        
       
    
if __name__ == "__main__":
    main()
    