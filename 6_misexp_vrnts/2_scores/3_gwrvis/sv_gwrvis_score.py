#!/bin/python
import argparse
import pandas as pd 
import numpy as np
from io import StringIO
from pybedtools import BedTool
from pathlib import Path

def main():
    # arguments passed 
    parser = argparse.ArgumentParser(description="")
    # required 
    parser.add_argument("--chrom", 
                        help ="chromosome", 
                        required = True)
    parser.add_argument("--gwrvis", 
                    help ="gwrvis scores", 
                    required = True)
    parser.add_argument("--sv_info", 
                help ="SV info", 
                required = True)
    parser.add_argument("--sv_bed", 
            help ="Bed file containing structural variants", 
            required = True)
    parser.add_argument("--out", 
        help ="Output directory", 
        required = True)
    
    global args
    args = parser.parse_args()
    chrom = args.chrom
    gwrvis_path = args.gwrvis
    sv_info_path = args.sv_info
    sv_bed_path = args.sv_bed
    out_dir = args.out
    
    out_dir_path = Path(out_dir)
    out_dir_path.mkdir(parents=True, exist_ok=True)

    chrom_num = int(chrom.split("chr")[1])

    data_type="gwrvis"
    print(f"Scoring SVs on {chrom} with {data_type} score ...")
    vrnts_bed_columns = {0:"chrom", 1:"start", 2:"end", 3:"vrnt_id"}

    bed_intersect_cols = {0:"vrnt_chrom", 1:"vrnt_start", 2:"vrnt_end", 3:"vrnt_id", 4:f"{data_type}_chrom", 
                          5:f"{data_type}_start", 6:f"{data_type}_end", 7:f"{data_type}_score", 8:"overlap"}

    # load SV info
    sv_info_df =pd.read_csv(sv_info_path, sep="\t", dtype={"plinkID": str}).rename(columns={"plinkID":"vrnt_id"})

    # load variants in windows 
    sv_bed = BedTool(sv_bed_path)
    sv_bed_sorted = sv_bed.sort()
    vrnts_in_windows_df = pd.read_csv(sv_bed_path, sep="\t", header=None, dtype={3:str}).rename(columns=vrnts_bed_columns)
    chrom_vrnts_in_windows_df = vrnts_in_windows_df[vrnts_in_windows_df.chrom == chrom_num]
    chrom_vrnts_in_windows = set(chrom_vrnts_in_windows_df.vrnt_id.unique())
    print(f"Number of variants in windows on {chrom}: {len(chrom_vrnts_in_windows)}")

    # load gwRVIS scores 
    gwrvis_bed = BedTool(gwrvis_path)

    # extract overlapping gwRVIS 
    sv_intersect_gwrvis_str = StringIO(str(sv_bed_sorted.intersect(gwrvis_bed, wo=True)))
    sv_intersect_gwrvis_df = pd.read_csv(sv_intersect_gwrvis_str, sep="\t", header=None, dtype={0: int, 3: str}).rename(columns=bed_intersect_cols)

    # only keep SV that intersect the correct chromosome 
    # gwRVIS bed files seem to have mix of chromosomes 
    sv_intersect_gwrvis_chrom_df = sv_intersect_gwrvis_df[sv_intersect_gwrvis_df.vrnt_chrom == chrom_num]
    # count variant with gwRVIS score 
    sv_with_gwrvis = set(sv_intersect_gwrvis_chrom_df.vrnt_id.unique())
    print(f"Number of variants with gwRVIS intersect {chrom}: {len(sv_with_gwrvis)}")

    if len(sv_with_gwrvis) > len(chrom_vrnts_in_windows):
        raise ValueError("Intersection file contains more variants than input file")
    if len(sv_with_gwrvis - chrom_vrnts_in_windows) != 0: 
        raise ValueError(f"Variants in intersection file not in input: {intersect_sv - chrom_vrnts_in_windows}")
    # expand dataframe so that each base pair in SV has a single score
    sv_intersect_gwrvis_chrom_expanded_df = sv_intersect_gwrvis_chrom_df.loc[sv_intersect_gwrvis_chrom_df.index.repeat(sv_intersect_gwrvis_chrom_df.overlap)][["vrnt_id", f"{data_type}_score"]]
    # check correct number of entries
    if sv_intersect_gwrvis_chrom_expanded_df.shape[0] != sv_intersect_gwrvis_chrom_df.overlap.sum():
        raise ValueError("Different number of entries and overlapping base pairs")

    vrnts_in_windows_gwrvis_df = chrom_vrnts_in_windows_df.copy()
    # calculate mean, meadian, first quartile, third quartile, max and  min gwRVIS per SV 
    gwrvis_metrics_df = pd.DataFrame(sv_intersect_gwrvis_chrom_expanded_df.groupby("vrnt_id")["gwrvis_score"].mean()).reset_index().rename(columns={"gwrvis_score":"mean"})
    gwrvis_metrics_df["median"] = sv_intersect_gwrvis_chrom_expanded_df.groupby("vrnt_id")["gwrvis_score"].median().tolist()
    gwrvis_metrics_df["quart1"] = sv_intersect_gwrvis_chrom_expanded_df.groupby("vrnt_id")["gwrvis_score"].quantile(0.25).tolist()
    gwrvis_metrics_df["quart3"] = sv_intersect_gwrvis_chrom_expanded_df.groupby("vrnt_id")["gwrvis_score"].quantile(0.75).tolist()
    # compute mean of 4 metrics 
    gwrvis_metrics_df["vitsios"] = gwrvis_metrics_df[["mean", "median", "quart1", "quart3"]].mean(axis=1)
    gwrvis_metrics_df["max"] = sv_intersect_gwrvis_chrom_expanded_df.groupby("vrnt_id")["gwrvis_score"].max().tolist()
    gwrvis_metrics_df["min"] = sv_intersect_gwrvis_chrom_expanded_df.groupby("vrnt_id")["gwrvis_score"].min().tolist()
    # add metrics
    vrnts_in_windows_gwrvis_df = pd.merge(vrnts_in_windows_gwrvis_df, 
                                          gwrvis_metrics_df, 
                                          on="vrnt_id", 
                                          how="left")
    # add structural variant types and length
    vrnts_in_windows_gwrvis_df = pd.merge(vrnts_in_windows_gwrvis_df, 
                                          sv_info_df[["vrnt_id", "SVLEN", "SVTYPE"]], 
                                          on="vrnt_id", 
                                          how="left")
    # check for NaNs
    print(f"NAs in dataframe: {vrnts_in_windows_gwrvis_df.isnull().values.any()}")
    print(f"Unique variants in final output: {len(vrnts_in_windows_gwrvis_df.vrnt_id.unique())}")
    # write results 
    path_out = out_dir_path.joinpath(f"{chrom}_sv_in_windows_{data_type}.tsv")
    vrnts_in_windows_gwrvis_df.to_csv(path_out, sep="\t", index=False) 

    
if __name__ == "__main__":
    main()