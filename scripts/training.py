#!/usr/bin/env python3
"""
Training script for learning statistical potentials for RNA 3D structures.

Implements Task 1:
 - all important constants are parameters (max distance, bin width, min sequence separation)
 - clean modular structure
 - command-line interface (--help)

Author: Negin Heidarifard
"""

import os
import argparse
import numpy as np
from Bio.PDB import PDBParser


BASES = {"A", "U", "C", "G"}


def parse_arguments():
    """Parse command-line parameters."""
    parser = argparse.ArgumentParser(
        description="Train statistical potentials from RNA PDB structures."
    )

    parser.add_argument(
        "--pdb_dir",
        required=True,
        help="Folder containing training PDB files (e.g. data/train_pdb)",
    )

    parser.add_argument(
        "--out_dir",
        required=True,
        help="Folder where score files will be written (e.g. scores)",
    )

    parser.add_argument(
        "--max_dist",
        type=float,
        default=20.0,
        help="Maximum distance cutoff in Å (default: 20.0)",
    )

    parser.add_argument(
        "--bin_width",
        type=float,
        default=1.0,
        help="Width of distance bins in Å (default: 1.0)",
    )

    parser.add_argument(
        "--min_sep",
        type=int,
        default=3,
        help="Minimum sequence separation (i -> i+min_sep); default: 3",
    )

    return parser.parse_args()


def get_bins(max_dist: float, bin_width: float) -> np.ndarray:
    """Return bin edges for distance histogram."""
    num_bins = int(max_dist / bin_width)
    return np.linspace(0.0, max_dist, num_bins + 1)


def train_potential(pdb_dir, out_dir, max_dist, bin_width, min_sep):
    """Main training routine."""

    parser = PDBParser(QUIET=True)
    bins = get_bins(max_dist, bin_width)
    num_bins = len(bins) - 1

    # 10 base-pair types
    pair_order = [
        "AA", "AU", "AC", "AG",
        "UU", "UC", "UG",
        "CC", "CG",
        "GG",
    ]
    pair_to_idx = {p: i for i, p in enumerate(pair_order)}

    # counts: 10 x num_bins
    obs_counts = np.zeros((len(pair_order), num_bins), dtype=float)
    # reference counts: XX pair
    ref_counts = np.zeros(num_bins, dtype=float)

    # --- loop over PDB files ---
    for fname in os.listdir(pdb_dir):
        if not fname.lower().endswith(".pdb"):
            continue

        path = os.path.join(pdb_dir, fname)
        print(f"Processing {path}")
        structure = parser.get_structure("RNA", path)
        model = structure[0]

        for chain in model:
            # collect residues in this chain that have C3' and are A/U/C/G
            residues = []
            for res in chain:
                if "C3'" not in res:
                    continue
                resname = res.get_resname().strip()[0]
                if resname not in BASES:
                    continue
                residues.append(res)

            n = len(residues)
            if n == 0:
                continue

            # pairwise distances with sequence separation >= min_sep
            for i in range(n):
                res_i = residues[i]
                base_i = res_i.get_resname().strip()[0]
                coord_i = res_i["C3'"].coord

                for j in range(i + min_sep, n):
                    res_j = residues[j]
                    base_j = res_j.get_resname().strip()[0]
                    coord_j = res_j["C3'"].coord

                    # Euclidean distance
                    d = np.linalg.norm(coord_i - coord_j)
                    if d > max_dist:
                        continue

                    # which bin?
                    bin_idx = int(d // bin_width)
                    if bin_idx < 0 or bin_idx >= num_bins:
                        continue

                    # base-pair type (sorted so AU == UA, etc.)
                    pair = "".join(sorted([base_i, base_j]))
                    if pair not in pair_to_idx:
                        continue

                    obs_counts[pair_to_idx[pair], bin_idx] += 1.0
                    ref_counts[bin_idx] += 1.0

    # --- convert counts to frequencies ---
    # avoid division by zero
    obs_sums = obs_counts.sum(axis=1, keepdims=True)
    obs_sums[obs_sums == 0.0] = 1.0
    obs_freq = obs_counts / obs_sums

    ref_sum = ref_counts.sum()
    if ref_sum == 0.0:
        ref_sum = 1.0
    ref_freq = ref_counts / ref_sum

    # --- compute pseudo-energy: -log(f_obs / f_ref) ---
    eps = 1e-8
    ratio = (obs_freq + eps) / (ref_freq + eps)
    scores = -np.log(ratio)

    # cap maximum value at +10 (as in TP)
    scores = np.minimum(scores, 10.0)

    # --- write score files ---
    os.makedirs(out_dir, exist_ok=True)
    for i, pair in enumerate(pair_order):
        out_path = os.path.join(out_dir, f"{pair}.txt")
        np.savetxt(out_path, scores[i], fmt="%.6f")
        print(f"  wrote {out_path}")

    print("Training completed. Score files written to:", out_dir)


if __name__ == "__main__":
    args = parse_arguments()
    train_potential(
        args.pdb_dir,
        args.out_dir,
        args.max_dist,
        args.bin_width,
        args.min_sep,
    )
