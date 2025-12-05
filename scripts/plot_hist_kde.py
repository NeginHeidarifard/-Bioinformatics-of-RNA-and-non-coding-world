#!/usr/bin/env python3
"""
Extra script: plot HISTOGRAM + KDE for one base-pair type.

- Recomputes C3'-C3' distances from training PDBs
- Builds a histogram with given bin width (or an automatic one)
- Overlays a simple Gaussian KDE

Useful to visualise:
 - effect of bin width
 - relation between histogram and KDE
"""

import os
import argparse
import numpy as np
from Bio.PDB import PDBParser
import matplotlib
matplotlib.use("Agg")  # WSL: no GUI
import matplotlib.pyplot as plt

BASES = {"A", "U", "C", "G"}


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot histogram + KDE of C3'-C3' distances for one base pair."
    )

    parser.add_argument(
        "--pdb_dir",
        required=True,
        help="Folder with training PDB files (e.g. data/train_pdb)",
    )

    parser.add_argument(
        "--out_path",
        required=True,
        help="Path of the output PNG (e.g. plots/hist_kde_AU.png)",
    )

    parser.add_argument(
        "--pair",
        default="XX",
        help="Base-pair type to analyse: AA, AU, ..., GG or XX for all bases pooled.",
    )

    parser.add_argument(
        "--max_dist",
        type=float,
        default=20.0,
        help="Maximum distance (Å) for the histogram (default: 20.0).",
    )

    parser.add_argument(
        "--bin_width",
        type=float,
        default=-1.0,
        help="Bin width (Å) of the histogram. If < 0, an automatic rule will be used.",
    )

    parser.add_argument(
        "--bandwidth",
        type=float,
        default=-1.0,
        help="Bandwidth (Å) for Gaussian KDE. If < 0, defaults to bin_width.",
    )

    return parser.parse_args()


def collect_distances(pdb_dir, max_dist, min_sep, pair):
    """
    Collect all C3'-C3' distances for the requested base pair.

    pair == "XX"  → ignore base identity (reference distribution)
    pair == "AU"  → only pairs whose sorted bases are "AU", etc.
    """
    parser = PDBParser(QUIET=True)
    all_d = []

    for fname in os.listdir(pdb_dir):
        if not fname.lower().endswith(".pdb"):
            continue

        path = os.path.join(pdb_dir, fname)
        structure = parser.get_structure("RNA", path)
        model = structure[0]

        for chain in model:
            residues = []
            for res in chain:
                if "C3'" not in res:
                    continue
                resname = res.get_resname().strip()[0]
                if resname not in BASES:
                    continue
                residues.append(res)

            n = len(residues)
            for i in range(n):
                res_i = residues[i]
                base_i = res_i.get_resname().strip()[0]
                coord_i = res_i["C3'"].coord

                for j in range(i + min_sep, n):
                    res_j = residues[j]
                    base_j = res_j.get_resname().strip()[0]
                    coord_j = res_j["C3'"].coord

                    d = np.linalg.norm(coord_i - coord_j)
                    if d > max_dist:
                        continue

                    if pair != "XX":
                        bp = "".join(sorted([base_i, base_j]))
                        if bp != pair:
                            continue

                    all_d.append(d)

    return np.array(all_d)


def gaussian_kde_1d(data, grid, bandwidth):
    """
    Simple 1D Gaussian KDE implemented with NumPy only.
    """
    if data.size == 0:
        return np.zeros_like(grid)

    diff = data[:, None] - grid[None, :]
    kernel = np.exp(-0.5 * (diff / bandwidth) ** 2)
    kde = kernel.mean(axis=0)
    # normalize so area ~1
    kde /= (kde.sum() * (grid[1] - grid[0]))
    return kde


def suggest_bin_widths(dists):
    """
    Suggest histogram bin widths using Scott's rule and Freedman–Diaconis rule.

    Returns (scott_width, fd_width).
    """
    n = dists.size
    if n == 0:
        return None, None

    std = np.std(dists)
    q25, q75 = np.percentile(dists, [25, 75])
    iqr = q75 - q25

    scott = 3.5 * std / (n ** (1.0 / 3.0)) if std > 0 else None
    fd = 2.0 * iqr / (n ** (1.0 / 3.0)) if iqr > 0 else None

    return scott, fd


def main():
    args = parse_arguments()

    max_dist = args.max_dist
    bin_width = args.bin_width
    bandwidth = args.bandwidth
    pair = args.pair.upper()

    print(f"Collecting distances for pair = {pair} ...")
    dists = collect_distances(
        args.pdb_dir,
        max_dist=max_dist,
        min_sep=3,        # same min_sep as in training
        pair=pair,
    )

    if dists.size == 0:
        print("No distances found for this pair; nothing to plot.")
        return

    print(f"Collected {dists.size} distances.")

    # --- choose bin_width automatically if requested ---
    if bin_width <= 0:
        scott, fd = suggest_bin_widths(dists)
        print("Suggested bin widths:")
        if scott is not None:
            print(f"  Scott's rule        ≈ {scott:.3f} Å")
        else:
            print("  Scott: n/a")
        if fd is not None:
            print(f"  Freedman–Diaconis   ≈ {fd:.3f} Å")
        else:
            print("  FD: n/a")

        if fd is not None and fd > 0:
            bin_width = fd
            print(f"Using FD bin width: {bin_width:.3f} Å")
        elif scott is not None and scott > 0:
            bin_width = scott
            print(f"Using Scott bin width: {bin_width:.3f} Å")
        else:
            bin_width = 1.0
            print("Fallback: using bin width = 1.0 Å")

    # --- choose bandwidth automatically if requested ---
    if bandwidth <= 0:
        bandwidth = bin_width
        print(f"Using bandwidth = bin_width = {bandwidth:.3f} Å")

    # --- build histogram ---
    bins = np.arange(0.0, max_dist + bin_width, bin_width)
    hist_counts, bin_edges = np.histogram(dists, bins=bins, density=True)
    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # --- build KDE on a fine grid ---
    grid = np.linspace(0.0, max_dist, 400)
    kde_vals = gaussian_kde_1d(dists, grid, bandwidth)

    # --- plot ---
    plt.figure(figsize=(6, 4))
    plt.bar(
        bin_centres,
        hist_counts,
        width=bin_width,
        alpha=0.5,
        label=f"Histogram (bin width = {bin_width:.2f} Å)",
        align="center",
    )

    plt.plot(
        grid,
        kde_vals,
        label=f"KDE (bandwidth = {bandwidth:.2f} Å)",
    )

    plt.xlabel("C3'-C3' distance (Å)")
    plt.ylabel("Probability density")
    title_pair = pair if pair != "XX" else "XX (all bases)"
    plt.title(f"Histogram + KDE for {title_pair}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out_path, dpi=150)
    plt.close()

    print("Wrote:", args.out_path)


if __name__ == "__main__":
    main()
