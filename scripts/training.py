#!/usr/bin/env python3
"""
Training script for learning statistical potentials for RNA 3D structures.

Implements:
- Task 1: constants as parameters (CLI, no hard-coded values)
- Task 4: alternative non-log scoring formulation (Postic et al., 2020)

Author: Negin Heidarifard
"""

import os
import argparse
import numpy as np
from Bio.PDB import PDBParser

BASES = {"A", "U", "C", "G"}


# -----------------------------
# Argument parsing
# -----------------------------
def parse_arguments():
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
        help="Folder where score files will be written",
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
        help="Minimum sequence separation (default: 3)",
    )

    parser.add_argument(
        "--formula",
        choices=["log", "nonlog"],
        default="log",
        help="Scoring formula: log (Sippl) or nonlog (information gain)",
    )

    parser.add_argument(
        "--epsilon",
        type=float,
        default=1e-8,
        help="Pseudocount to avoid division by zero (default: 1e-8)",
    )

    parser.add_argument(
        "--cap",
        type=float,
        default=10.0,
        help="Absolute cap for pseudo-energy values (default: 10.0)",
    )

    return parser.parse_args()


# -----------------------------
# Utility functions
# -----------------------------
def get_bins(max_dist: float, bin_width: float) -> np.ndarray:
    """Return bin edges for distance histogram."""
    num_bins = int(max_dist / bin_width)
    return np.linspace(0.0, max_dist, num_bins + 1)


# -----------------------------
# Main training routine
# -----------------------------
def train_potential(args):
    parser = PDBParser(QUIET=True)
    bins = get_bins(args.max_dist, args.bin_width)
    num_bins = len(bins) - 1

    # 10 unordered nucleotide pair types
    pair_order = [
        "AA", "AU", "AC", "AG",
        "UU", "UC", "UG",
        "CC", "CG",
        "GG",
    ]
    pair_to_idx = {p: i for i, p in enumerate(pair_order)}

    obs_counts = np.zeros((len(pair_order), num_bins), dtype=float)
    ref_counts = np.zeros(num_bins, dtype=float)

    # -----------------------------
    # Loop over PDB files
    # -----------------------------
    for fname in os.listdir(args.pdb_dir):
        if not fname.lower().endswith(".pdb"):
            continue

        path = os.path.join(args.pdb_dir, fname)
        print(f"Processing {path}")

        structure = parser.get_structure("RNA", path)
        model = structure[0]

        for chain in model:
            residues = []
            for res in chain:
                if "C3'" not in res:
                    continue
                base = res.get_resname().strip()[0]
                if base not in BASES:
                    continue
                residues.append(res)

            n = len(residues)
            if n == 0:
                continue

            for i in range(n):
                res_i = residues[i]
                base_i = res_i.get_resname().strip()[0]
                coord_i = res_i["C3'"].coord

                for j in range(i + args.min_sep, n):
                    res_j = residues[j]
                    base_j = res_j.get_resname().strip()[0]
                    coord_j = res_j["C3'"].coord

                    d = np.linalg.norm(coord_i - coord_j)
                    if d > args.max_dist:
                        continue

                    bin_idx = int(d // args.bin_width)
                    if bin_idx < 0 or bin_idx >= num_bins:
                        continue

                    pair = "".join(sorted([base_i, base_j]))
                    if pair not in pair_to_idx:
                        continue

                    obs_counts[pair_to_idx[pair], bin_idx] += 1.0
                    ref_counts[bin_idx] += 1.0

    # -----------------------------
    # Frequencies
    # -----------------------------
    obs_sums = obs_counts.sum(axis=1, keepdims=True)
    obs_sums[obs_sums == 0.0] = 1.0
    obs_freq = obs_counts / obs_sums

    ref_sum = ref_counts.sum()
    if ref_sum == 0.0:
        ref_sum = 1.0
    ref_freq = ref_counts / ref_sum

    # -----------------------------
    # Scoring
    # -----------------------------
    eps = args.epsilon

    if args.formula == "log":
        scores = -np.log((obs_freq + eps) / (ref_freq + eps))
    elif args.formula == "nonlog":
        scores = -((obs_freq - ref_freq) / (ref_freq + eps))
    else:
        raise ValueError("Unknown scoring formula")

    scores = np.clip(scores, -args.cap, args.cap)

    # -----------------------------
    # Output
    # -----------------------------
    os.makedirs(args.out_dir, exist_ok=True)

    for i, pair in enumerate(pair_order):
        out_path = os.path.join(args.out_dir, f"{pair}.txt")
        np.savetxt(out_path, scores[i], fmt="%.6f")
        print(f"  wrote {out_path}")

    print("Training completed.")
    print(f"Formula: {args.formula}")
    print("Output directory:", args.out_dir)


# -----------------------------
# Entry point
# -----------------------------
if __name__ == "__main__":
    args = parse_arguments()
    train_potential(args)
