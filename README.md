````markdown
Statistical scoring of RNA 3D structures
This project is part of the RNA structure bioinformatics practical course (Master 2 GENIOMHE, Université Paris-Saclay). It implements a simple statistical potential to approximate the Gibbs free energy of RNA 3D conformations using inter-atomic distance distributions between ribonucleotides.

## Project goals

The aim is to build an “objective function” that assigns a pseudo-energy to an RNA structure.
The function is trained on experimentally determined RNA structures and is later used to score new conformations, for example models from the RNA-Puzzles dataset.

Concretely, the project:

- computes distance distributions between C3' atoms in native RNA structures;
- derives a pairwise statistical potential for the ten nucleotide pairs (AA, AU, AC, AG, UU, UC, UG, CC, CG, GG);
- uses this potential to estimate the Gibbs free energy of a candidate RNA conformation.

## Repository structure

```text
.
├── data/
│   ├── train_pdb/          # training RNA PDB files (native structures)
│   └── test_pdb/           # RNA-Puzzles or other test structures
├── scores/
│   └── *.txt               # learned scoring profiles for each base pair
├── scripts/
│   ├── training.py
│   ├── plot_profiles.py    # or plot_profiles.R
│   └── score_structure.py
└── README.md
````

The `data/` folders are not tracked here by default; local paths can be adapted in the command-line arguments of the scripts.

## Method overview

### 1. Training the statistical potential

The script `training.py` learns pairwise distance distributions from a dataset of RNA structures.

**Atom selection**

* Only C3' atoms of ribonucleotides are considered.
* Only intra-chain contacts are used.
* Residue pairs closer than three positions along the sequence (i, i+1, i+2, i+3) are discarded to remove trivial local geometry.

**Distance binning**

* For each residue pair type (AA, AU, AC, AG, UU, UC, UG, CC, CG, GG), all C3'–C3' distances are computed.
* Distances are placed into 20 bins between 0 and 20 Å (bin width 1 Å).

**Observed and reference probabilities**

* Observed probability `f_obs(i,j,r)` is computed by normalising counts over all bins.
* Reference distribution `f_ref(r)` is computed by pooling all bases as “X–X”.

**Score (pseudo-energy) calculation**

* Score is defined as:

  `s(i,j,r) = - log( f_obs(i,j,r) / f_ref(r) )`

* Scores are capped at +10 (TP instructions).

**Output**

* Ten text files are produced, one per nucleotide pair.
* Each with 20 scoring values (one per distance bin).

---

### 2. Plotting interaction profiles

The script `plot_profiles.py` visualises the learned potentials:

* Reads score files from `scores/`
* Plots energy vs distance bin centre (0.5, 1.5, …)
* Saves PNG or PDF figures

---

### 3. Scoring new RNA structures

The script `score_structure.py`:

1. Parses PDB and extracts C3' atoms.
2. Computes all valid intra-chain distances.
3. Determines base-pair type.
4. Uses linear interpolation between bins.
5. Sums all contributions to estimate Gibbs free energy.

Lower total energy = more native-like conformation.

---

## Command-line usage (planned)

```bash
# Train the potential
python scripts/training.py \
    --pdb_dir data/train_pdb \
    --max_dist 20.0 \
    --bin_width 1.0 \
    --out_dir scores/

# Plot interaction profiles
python scripts/plot_profiles.py \
    --scores_dir scores/ \
    --out_dir plots/

# Score a candidate structure
python scripts/score_structure.py \
    --pdb data/test_pdb/model.pdb \
    --scores_dir scores/
```

## Implementation notes

* Built using Python + Biopython for PDB parsing.
* Code modularised into small reusable functions.
* Each script will include a `--help` option.
* Git is used for version control and development tracking.
