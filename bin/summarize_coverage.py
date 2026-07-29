#!/usr/bin/env python3
import sys
import os
import argparse

__version__ = "1.0.0"

def summarize_depth(depth_file, sample_id, windowsize=500):
    depths = []
    if os.path.exists(depth_file):
        with open(depth_file, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    try:
                        depths.append(int(parts[2]))
                    except ValueError:
                        continue

    n = len(depths)
    if n == 0:
        return {
            "sample": sample_id,
            "mean_depth": 0.0,
            "min_depth": 0,
            "max_depth": 0,
            "percent_above_100x": 0.0,
            "percent_above_1000x": 0.0,
            "percent_above_5000x": 0.0
        }

    sum_depth = sum(depths)
    mean_depth = sum_depth / n
    min_depth = min(depths)
    max_depth = max(depths)
    above100 = sum(1 for d in depths if d >= 100)
    above1000 = sum(1 for d in depths if d >= 1000)
    above5000 = sum(1 for d in depths if d >= 5000)

    return {
        "sample": sample_id,
        "mean_depth": round(mean_depth, 1),
        "min_depth": min_depth,
        "max_depth": max_depth,
        "percent_above_100x": round((above100 / n) * 100, 1),
        "percent_above_1000x": round((above1000 / n) * 100, 1),
        "percent_above_5000x": round((above5000 / n) * 100, 1)
    }

def main():
    parser = argparse.ArgumentParser(description="Summarize samtools depth output.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--depth-file", required=True, help="Input samtools depth TSV")
    parser.add_argument("--sample-id", required=True, help="Sample identifier")
    parser.add_argument("--output-tsv", required=True, help="Output summary TSV file")
    parser.add_argument("--windowsize", type=int, default=500, help="Window size for rolling average")

    args = parser.parse_args()

    stats = summarize_depth(args.depth_file, args.sample_id, args.windowsize)

    os.makedirs(os.path.dirname(os.path.abspath(args.output_tsv)), exist_ok=True)
    with open(args.output_tsv, 'w', encoding='utf-8') as f:
        f.write("sample\tmean_depth\tmin_depth\tmax_depth\tpercent_above_100x\tpercent_above_1000x\tpercent_above_5000x\n")
        f.write(f"{stats['sample']}\t{stats['mean_depth']}\t{stats['min_depth']}\t{stats['max_depth']}\t{stats['percent_above_100x']}\t{stats['percent_above_1000x']}\t{stats['percent_above_5000x']}\n")

    print(f"Coverage summary written to {args.output_tsv}")

if __name__ == "__main__":
    main()
