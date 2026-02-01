#!/usr/bin/env python3
"""
Score an RNA 3D structure using an all-atom (base-only) statistical potential.

This script is the all-atom counterpart of score_structure.py and is consistent
with training_all_atom.py.

Given:
 - a PDB file (RNA structure),
 - a folder of all-atom score profiles (one file per atom-pair type),

the script:
 - extracts base atoms (excluding backbone atoms),
 - computes all valid atom–atom distances (intrachain, residue separation >= min_sep),
 - retrieves pseudo-energy values by linear interpolation,
 - sums all contributions to obtain an estimated pseudo–Gibbs free energy.

Author: Negin Heidarifard
"""

import os
import argparse
import numpy as np
from Bio.PDB import PDBParser


BASES = {"A", "U", "C", "G"}

# Backbone atoms to exclude (same as in training_all_atom.py)
BACKBONE_ATOMS = {
    "P", "OP1", "OP2",
    "O5'", "C5'", "C4'", "O4'",
    "C3'", "O3'", "C2'", "C1'",
    "O2'"
}


# ---------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Score an RNA structure using an all-atom statistical potential."
    )

    parser.add_argument(
        "--pdb",
        required=True,
        help="Path to the RNA PDB file to score",
    )

    parser.add_argument(
        "--scores_dir",
        required=True,
        help="Folder containing all-atom score profiles (*.txt)",
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
        help="Minimum residue separation (i -> i+min_sep); default: 3",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------
def residue_base(res):
    resname = res.get_resname().strip()
    if resname and resname[0] in BASES:
        return resname[0]
    return None


def extract_base_atoms(res):
    """
    Return a list of (atom_type, coord) for base atoms in a residue.
    atom_type is of the form <Base>_<AtomName>, e.g. A_N1.
    """
    base = residue_base(res)
    if base is None:
        return []

    atoms = []
    for atom in res:
        name = atom.get_name().strip()
        if name in BACKBONE_ATOMS:
            continue
        atom_type = f"{base}_{name}"
        atoms.append((atom_type, atom.coord))

    return atoms


def pair_key(t1, t2):
    """Canonical atom-pair key (order-independent)."""
    return "__".join(sorted([t1, t2]))


def load_score_profiles(scores_dir):
    """
    Load all score profiles from a directory.
    Each file must be named <atom1>__<atom2>.txt
    """
    profiles = {}

    for fname in os.listdir(scores_dir):
        if not fname.endswith(".txt"):
            continue
        key = fname.replace(".txt", "")
        path = os.path.join(scores_dir, fname)
        profiles[key] = np.loadtxt(path)

    if not profiles:
        raise RuntimeError(f"No score files found in {scores_dir}")

    return profiles


def interpolate_score(scores_array, d, bin_width, max_dist):
    """
    Linearly interpolate the pseudo-energy for distance d.
    """
    if d < 0.0 or d > max_dist:
        return None

    num_bins = scores_array.shape[0]
    idx = d / bin_width

    low = int(np.floor(idx))
    high = low + 1

    if low < 0:
        return scores_array[0]
    if high >= num_bins:
        return scores_array[-1]

    t = idx - low
    return (1.0 - t) * scores_array[low] + t * scores_array[high]


# ---------------------------------------------------------------------
# Main scoring routine
# ---------------------------------------------------------------------
def score_structure(pdb_path, scores_dir, max_dist, bin_width, min_sep):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", pdb_path)
    model = structure[0]

    profiles = load_score_profiles(scores_dir)

    total_score = 0.0
    num_pairs = 0

    for chain in model:
        # collect residues with valid bases
        residues = [r for r in chain if residue_base(r)]
        atoms_per_residue = [extract_base_atoms(r) for r in residues]

        n = len(residues)
        if n == 0:
            continue

        for i in range(n):
            atoms_i = atoms_per_residue[i]
            if not atoms_i:
                continue

            for j in range(i + min_sep, n):
                atoms_j = atoms_per_residue[j]
                if not atoms_j:
                    continue

                for (t1, c1) in atoms_i:
                    for (t2, c2) in atoms_j:
                        d = np.linalg.norm(c1 - c2)
                        if d > max_dist:
                            continue

                        key = pair_key(t1, t2)
                        if key not in profiles:
                            continue

                        s = interpolate_score(
                            profiles[key], d, bin_width, max_dist
                        )
                        if s is None:
                            continue

                        total_score += s
                        num_pairs += 1

    return total_score, num_pairs


# ---------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------
def main():
    args = parse_arguments()

    print(f"Scoring structure: {args.pdb}")
    print(f"Using all-atom profiles in: {args.scores_dir}")
    print(
        f"max_dist = {args.max_dist} Å, "
        f"bin_width = {args.bin_width} Å, "
        f"min_sep = {args.min_sep}"
    )

    total_score, num_pairs = score_structure(
        args.pdb,
        args.scores_dir,
        args.max_dist,
        args.bin_width,
        args.min_sep,
    )

    if num_pairs == 0:
        print("Warning: no valid atom–atom pairs found.")
    else:
        print(f"Number of contributing atom pairs: {num_pairs}")
        print(
            f"Estimated pseudo-Gibbs free energy: "
            f"{total_score:.4f} (arbitrary units)"
        )


if __name__ == "__main__":
    main()
