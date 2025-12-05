"""
Third script: score an RNA 3D structure with the learned statistical potential.

Given:
 - a PDB file (RNA structure),
 - a folder of score profiles produced by training.py,

the script:
 - computes all valid C3'-C3' distances (intrachain, seq. separation >= min_sep),
 - for each pair, retrieves a pseudo-energy value by linear interpolation
   in the corresponding base-pair profile,
 - sums all contributions to obtain an estimated Gibbs free energy.

Author: Negin Heidarifard
"""

import os
import argparse
import numpy as np
from Bio.PDB import PDBParser


BASES = {"A", "U", "C", "G"}

PAIR_ORDER = [
    "AA", "AU", "AC", "AG",
    "UU", "UC", "UG",
    "CC", "CG",
    "GG",
]
PAIR_TO_IDX = {p: i for i, p in enumerate(PAIR_ORDER)}


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Score an RNA 3D structure using the learned statistical potential."
    )

    parser.add_argument(
        "--pdb",
        required=True,
        help="Path to the RNA PDB file to score (e.g. data/test_pdb/model.pdb)",
    )

    parser.add_argument(
        "--scores_dir",
        required=True,
        help="Folder containing score profiles (e.g. scores)",
    )

    parser.add_argument(
        "--max_dist",
        type=float,
        default=20.0,
        help="Maximum distance cutoff in Å (must match training; default: 20.0)",
    )

    parser.add_argument(
        "--bin_width",
        type=float,
        default=1.0,
        help="Bin width in Å used during training (default: 1.0)",
    )

    parser.add_argument(
        "--min_sep",
        type=int,
        default=3,
        help="Minimum sequence separation (i -> i+min_sep); default: 3",
    )

    return parser.parse_args()


def load_score_profiles(scores_dir):
    """
    Load score profiles for all base-pair types into a dictionary:
    pair -> 1D numpy array of length num_bins.
    """
    profiles = {}
    for pair in PAIR_ORDER:
        path = os.path.join(scores_dir, f"{pair}.txt")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Score file not found for {pair}: {path}")
        profiles[pair] = np.loadtxt(path)
    return profiles


def interpolate_score(scores_array, d, bin_width, max_dist):
    """
    Linearly interpolate the pseudo-energy for a given distance d
    using a 1D array of scores defined on contiguous bins of width bin_width
    from 0 to max_dist.

    Bins are assumed to cover [0, bin_width), [bin_width, 2*bin_width), ..., [max_dist-bin_width, max_dist).

    We convert d to a fractional index and linearly interpolate between
    scores[low] and scores[high].
    """
    if d < 0.0 or d > max_dist:
        return None  # caller should ignore

    num_bins = scores_array.shape[0]

    # fractional index
    idx = d / bin_width
    low = int(np.floor(idx))
    high = low + 1

    if low < 0:
        return scores_array[0]
    if high >= num_bins:
        return scores_array[-1]

    t = idx - low  # between 0 and 1
    return (1.0 * (1.0 - t)) * scores_array[low] + (1.0 * t) * scores_array[high]


def score_structure(pdb_path, scores_dir, max_dist, bin_width, min_sep):
    """
    Main scoring routine: returns (total_score, num_pairs).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", pdb_path)
    model = structure[0]

    profiles = load_score_profiles(scores_dir)

    total_score = 0.0
    num_pairs = 0

    for chain in model:
        # collect residues with C3' and standard bases
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

                pair = "".join(sorted([base_i, base_j]))
                if pair not in profiles:
                    continue

                s = interpolate_score(profiles[pair], d, bin_width, max_dist)
                if s is None:
                    continue

                total_score += s
                num_pairs += 1

    return total_score, num_pairs


def main():
    args = parse_arguments()

    print(f"Scoring structure: {args.pdb}")
    print(f"Using score profiles in: {args.scores_dir}")
    print(f"max_dist = {args.max_dist} Å, bin_width = {args.bin_width} Å, min_sep = {args.min_sep}")

    total_score, num_pairs = score_structure(
        args.pdb,
        args.scores_dir,
        args.max_dist,
        args.bin_width,
        args.min_sep,
    )

    if num_pairs == 0:
        print("Warning: no valid C3'-C3' pairs found under the given constraints.")
    print(f"Number of contributing pairs: {num_pairs}")
    print(f"Estimated pseudo-Gibbs free energy: {total_score:.4f} (arbitrary units)")


if __name__ == "__main__":
    main()
