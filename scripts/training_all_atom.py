#!/usr/bin/env python3
"""
All-atom (base-only) training script for RNA statistical potentials.

Learns distance-dependent potentials for atom-atom pairs:
u_ab(r) = -log( f_obs_ab(r) / f_ref_XX(r) )   (or nonlog variant)

- intrachain only
- min sequence separation >= min_sep
- max distance cutoff
- bin width
- base atoms only (backbone excluded)

Author: Negin Heidarifard
"""

import os
import argparse
from collections import defaultdict
import numpy as np
from Bio.PDB import PDBParser

BASES = {"A", "U", "C", "G"}

BACKBONE_ATOMS = {
    "P", "OP1", "OP2",
    "O5'", "C5'", "C4'", "O4'",
    "C3'", "O3'", "C2'", "C1'",
    "O2'"
}


def parse_args():
    p = argparse.ArgumentParser(description="Train all-atom (base-only) statistical potentials from RNA PDBs.")
    p.add_argument("--pdb_dir", required=True, help="Folder containing training PDB files.")
    p.add_argument("--out_dir", required=True, help="Output folder for learned score profiles.")
    p.add_argument("--max_dist", type=float, default=20.0, help="Max distance cutoff in Å (default: 20).")
    p.add_argument("--bin_width", type=float, default=1.0, help="Bin width in Å (default: 1).")
    p.add_argument("--min_sep", type=int, default=3, help="Minimum sequence separation (default: 3).")
    p.add_argument("--formula", choices=["log", "nonlog"], default="log",
                   help="Scoring formula: log (Sippl) or nonlog (information gain).")
    p.add_argument("--epsilon", type=float, default=1e-8, help="Pseudocount to avoid div by zero.")
    p.add_argument("--cap", type=float, default=10.0, help="Cap for pseudo-energy values (default: 10).")
    return p.parse_args()


def get_bin_edges(max_dist: float, bin_width: float) -> np.ndarray:
    n_bins = int(max_dist / bin_width)
    return np.linspace(0.0, max_dist, n_bins + 1)


def residue_base(res) -> str | None:
    rn = res.get_resname().strip()
    if not rn:
        return None
    b = rn[0]
    return b if b in BASES else None


def extract_base_atoms(res):
    """
    Return list of (atom_type, coord) for base atoms in residue.
    atom_type = <Base>_<AtomName>
    """
    b = residue_base(res)
    if b is None:
        return []

    atoms = []
    for atom in res:
        name = atom.get_name().strip()
        if name in BACKBONE_ATOMS:
            continue
        atoms.append((f"{b}_{name}", atom.coord))
    return atoms


def pair_key(t1: str, t2: str) -> str:
    """Canonical unordered key for an atom-type pair."""
    return "__".join(sorted([t1, t2]))


def compute_scores(obs_counts: np.ndarray, ref_counts: np.ndarray, epsilon: float, cap: float, formula: str) -> np.ndarray:
    """
    obs_counts: (n_pairs, n_bins)
    ref_counts: (n_bins,)
    returns: scores (n_pairs, n_bins)
    """
    # frequencies
    obs_sums = obs_counts.sum(axis=1, keepdims=True)
    obs_sums[obs_sums == 0.0] = 1.0
    f_obs = obs_counts / obs_sums

    ref_sum = ref_counts.sum()
    if ref_sum == 0.0:
        ref_sum = 1.0
    f_ref = ref_counts / ref_sum

    ratio = (f_obs + epsilon) / (f_ref + epsilon)

    if formula == "log":
        scores = -np.log(ratio)
    else:
        # "nonlog" option: information gain-like linear version
        # (kept simple and consistent with your earlier extension)
        scores = 1.0 - ratio

    # cap upper bound (as TP style)
    scores = np.minimum(scores, cap)
    return scores


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    bins = get_bin_edges(args.max_dist, args.bin_width)
    n_bins = len(bins) - 1

    parser = PDBParser(QUIET=True)

    # We don't know all atom-pair types ahead of time -> use dict then freeze to arrays
    pair_to_counts = defaultdict(lambda: np.zeros(n_bins, dtype=float))
    ref_counts = np.zeros(n_bins, dtype=float)

    pdb_files = [f for f in os.listdir(args.pdb_dir) if f.lower().endswith(".pdb")]
    pdb_files.sort()

    print("Training all-atom potential")
    print("pdb_dir:", args.pdb_dir)
    print("out_dir:", args.out_dir)
    print(f"max_dist={args.max_dist}, bin_width={args.bin_width}, min_sep={args.min_sep}, formula={args.formula}")
    print("Number of PDB files:", len(pdb_files))

    for fname in pdb_files:
        path = os.path.join(args.pdb_dir, fname)
        print("Processing", path)

        structure = parser.get_structure("RNA", path)
        model = structure[0]

        for chain in model:
            # Collect residues (keep ordering) and their base atoms
            residues = []
            atoms_per_res = []
            for res in chain:
                if residue_base(res) is None:
                    continue
                base_atoms = extract_base_atoms(res)
                if not base_atoms:
                    continue
                residues.append(res)
                atoms_per_res.append(base_atoms)

            n = len(residues)
            if n == 0:
                continue

            # Pair residues with sequence sep >= min_sep
            for i in range(n):
                atoms_i = atoms_per_res[i]
                for j in range(i + args.min_sep, n):
                    atoms_j = atoms_per_res[j]

                    # Atom-atom distances between the two residues
                    for (t1, c1) in atoms_i:
                        for (t2, c2) in atoms_j:
                            d = np.linalg.norm(c1 - c2)
                            if d > args.max_dist:
                                continue

                            bin_idx = int(d // args.bin_width)
                            if bin_idx < 0 or bin_idx >= n_bins:
                                continue

                            k = pair_key(t1, t2)
                            pair_to_counts[k][bin_idx] += 1.0
                            ref_counts[bin_idx] += 1.0

    # Freeze dict to arrays for score computation
    pair_keys = sorted(pair_to_counts.keys())
    obs_counts = np.stack([pair_to_counts[k] for k in pair_keys], axis=0) if pair_keys else np.zeros((0, n_bins))

    scores = compute_scores(obs_counts, ref_counts, args.epsilon, args.cap, args.formula)

    # Save reference distribution too (useful for debugging / reproducibility)
    ref_path = os.path.join(args.out_dir, "XX.txt")
    np.savetxt(ref_path, ref_counts / max(ref_counts.sum(), 1.0), fmt="%.8f")
    print("Wrote reference distribution:", ref_path)

    # Save per atom-pair score profiles
    for idx, k in enumerate(pair_keys):
        out_path = os.path.join(args.out_dir, f"{k}.txt")
        np.savetxt(out_path, scores[idx], fmt="%.6f")

    print("Done.")
    print("Number of atom-pair types learned:", len(pair_keys))
    print("Example files:", (pair_keys[:5] if len(pair_keys) >= 5 else pair_keys))


if __name__ == "__main__":
    main()
