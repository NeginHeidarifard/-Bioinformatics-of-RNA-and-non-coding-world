#!/usr/bin/env python3
"""
Second script: plot interaction profiles (score vs distance).

Reads the score files produced by training.py (one file per base pair),
reconstructs the distance axis from max_dist and bin_width, and saves
one plot per base pair in the output folder.
"""

import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")  # WSL: no GUI
import matplotlib.pyplot as plt

PAIR_ORDER = [
    "AA", "AU", "AC", "AG",
    "UU", "UC", "UG",
    "CC", "CG",
    "GG",
]


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot statistical potential profiles for RNA base pairs."
    )

    parser.add_argument(
        "--scores_dir",
        required=True,
        help="Folder containing score files (e.g. scores)",
    )

    parser.add_argument(
        "--out_dir",
        required=True,
        help="Folder where plots will be saved (e.g. plots)",
    )

    parser.add_argument(
        "--max_dist",
        type=float,
        default=20.0,
        help="Maximum distance used during training (Å).",
    )

    parser.add_argument(
        "--bin_width",
        type=float,
        default=1.0,
        help="Bin width used during training (Å).",
    )

    return parser.parse_args()


def plot_profiles(scores_dir, out_dir, max_dist, bin_width):
    os.makedirs(out_dir, exist_ok=True)

    # distance axis: bin centres
    num_bins = int(max_dist / bin_width)
    bin_centres = np.arange(num_bins) * bin_width + (bin_width / 2.0)

    for pair in PAIR_ORDER:
        path = os.path.join(scores_dir, f"{pair}.txt")
        if not os.path.exists(path):
            print(f"[warning] score file not found for {pair}: {path}")
            continue

        scores = np.loadtxt(path)

        if scores.size == 0:
            print(f"[warning] empty score file for {pair}")
            continue

        # make sure length matches num_bins
        if scores.shape[0] != num_bins:
            print(f"[warning] {pair}: expected {num_bins} bins, "
                  f"found {scores.shape[0]}")
            local_centres = (
                np.arange(scores.shape[0]) * bin_width + (bin_width / 2.0)
            )
        else:
            local_centres = bin_centres

        plt.figure()
        plt.plot(local_centres, scores, marker="o", linewidth=2, markersize=4)
        plt.xlabel("C3'-C3' distance (Å)", fontsize=11)
        plt.ylabel("Pseudo-energy (a.u.)", fontsize=11)
        plt.title(f"Interaction profile: {pair}", fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        out_path = os.path.join(out_dir, f"{pair}.png")
        plt.savefig(out_path, dpi=150)
        plt.close()
        print(f"  wrote {out_path}")

    print("All profiles plotted.")


if __name__ == "__main__":
    args = parse_arguments()
    plot_profiles(
        args.scores_dir,
        args.out_dir,
        args.max_dist,
        args.bin_width,
    )
