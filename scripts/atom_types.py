#!/usr/bin/env python3
"""
Extract base atom types from RNA PDB structures (All-atom Task 2, Step 2.1).

Atom types are defined as: <Base>_<AtomName>
Example: A_N1, G_O6, U_N3

Backbone atoms are explicitly excluded.

Author: Negin Heidarifard
"""

import os
from Bio.PDB import PDBParser

BASES = {"A", "U", "C", "G"}

BACKBONE_ATOMS = {
    "P", "OP1", "OP2",
    "O5'", "C5'", "C4'", "O4'",
    "C3'", "O3'", "C2'", "C1'",
    "O2'"   
    
}



def extract_atom_types(pdb_dir):
    parser = PDBParser(QUIET=True)
    atom_types = set()

    for fname in os.listdir(pdb_dir):
        if not fname.lower().endswith(".pdb"):
            continue

        path = os.path.join(pdb_dir, fname)
        structure = parser.get_structure("RNA", path)
        model = structure[0]

        for chain in model:
            for res in chain:
                resname = res.get_resname().strip()
                if not resname or resname[0] not in BASES:
                    continue

                base = resname[0]

                for atom in res:
                    atom_name = atom.get_name().strip()
                    if atom_name in BACKBONE_ATOMS:
                        continue

                    atom_type = f"{base}_{atom_name}"
                    atom_types.add(atom_type)

    return sorted(atom_types)


if __name__ == "__main__":
    pdb_dir = "data/train_pdb"
    atom_types = extract_atom_types(pdb_dir)

    print("Number of base atom types:", len(atom_types))
    print("Atom types:")
    for at in atom_types:
        print(at)


