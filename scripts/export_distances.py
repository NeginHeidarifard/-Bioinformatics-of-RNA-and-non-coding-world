#!/usr/bin/env python3
"""
Export C3'-C3' distance distributions to CSV for KDE analysis (Task 3).

This script extracts distances from training RNA PDB files and writes them
to a CSV file for further analysis in R.

Author: Negin Heidarifard
"""

import os
import argparse
import csv
import numpy as np
from Bio.PDB import PDBParser

BASES = {"A", "U", "C", "G"}


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Export C3'-C3' distances to CSV for KDE analysis"
    )
    parser.add_argument("--pdb_dir", required=True, help="Training PDB directory")
    parser.add_argument("--pair", default="AU",
                        help="Base pair to extract (e.g. AU or XX for all)")
    parser.add_argument("--max_dist", type=float, default=20.0,
                        help="Maximum distance cutoff (Å)")
    parser.add_argument("--min_sep", type=int, default=3,
                        help="Minimum sequence separation")
    parser.add_argument("--out_csv", required=True,
                        help="Output CSV file")
    return parser.parse_args()


def main():
    args = parse_arguments()
    parser = PDBParser(QUIET=True)

    distances = []

    for fname in os.listdir(args.pdb_dir):
        if not fname.lower().endswith(".pdb"):
            continue

        path = os.path.join(args.pdb_dir, fname)
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
            for i in range(n):
                ri = residues[i]
                bi = ri.get_resname().strip()[0]
                ci = ri["C3'"].coord

                for j in range(i + args.min_sep, n):
                    rj = residues[j]
                    bj = rj.get_resname().strip()[0]
                    cj = rj["C3'"].coord

                    d = np.linalg.norm(ci - cj)
                    if d > args.max_dist:
                        continue

                    pair = "".join(sorted([bi, bj]))
                    if args.pair != "XX" and pair != args.pair:
                        continue

                    distances.append(d)

    with open(args.out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["distance"])
        for d in distances:
            writer.writerow([d])

    print(f"Wrote {len(distances)} distances to {args.out_csv}")


if __name__ == "__main__":
    main()
