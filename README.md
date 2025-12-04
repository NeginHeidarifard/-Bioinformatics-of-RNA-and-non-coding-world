Statistical scoring of RNA 3D structures
This project is part of the RNA structure bioinformatics practical course (Master 2 GENIOMHE, Université Paris‑Saclay). It implements a simple statistical potential to approximate the Gibbs free energy of RNA 3D conformations using inter‑atomic distance distributions between ribonucleotides.​

Project goals
The aim is to build an “objective function” that assigns a pseudo‑energy to an RNA structure.​
The function is trained on experimentally determined RNA structures and is later used to score new conformations, for example models from the RNA‑Puzzles dataset.​

Concretely, the project:

computes distance distributions between C3' atoms in native RNA structures;​

derives a pairwise statistical potential for the ten nucleotide pairs (AA, AU, AC, AG, UU, UC, UG, CC, CG, GG);​

uses this potential to estimate the Gibbs free energy of a candidate RNA conformation.​

Repository structure
text
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
The data/ folders are not tracked here by default; local paths can be adapted in the command‑line arguments of the scripts.​

Method overview
1. Training the statistical potential
The script training.py learns pairwise distance distributions from a dataset of RNA structures.​

Atom selection

Only C3' atoms of ribonucleotides are considered.​

Only intra‑chain contacts are used.​

Residue pairs closer than three positions along the sequence (i, i+1, i+2, i+3) are discarded to remove trivial local geometry.​

Distance binning

For each residue pair type (AA, AU, AC, AG, UU, UC, UG, CC, CG, GG), all C3'–C3' distances are computed.​

Distances are placed into 20 bins between 0 and 20 Å (bin width 1 Å).​

Observed and reference probabilities

For each pair type and distance bin r, an observed frequency f_obs(i,j,r) is computed by normalising the counts over all bins.​

A reference distribution f_ref(r) is obtained in the same way but ignoring nucleotide identities (all bases pooled as “X–X”).​

Score (pseudo‑energy) calculation

For each nucleotide pair and bin, the score is defined as:
s(i,j,r) = - log( f_obs(i,j,r) / f_ref(r) )
which follows the definition of the pseudo‑energy in the TP.​

Scores are capped at +10, following the TP instructions.​

Output

Ten text files are produced, one per nucleotide pair.​

Each file contains 20 lines, corresponding to the 20 distance bins, with a single scoring value per line.​

2. Plotting interaction profiles
The script plot_profiles.py (or plot_profiles.R) visualises the learned potentials.​

It reads the score files for all nucleotide pairs from scores/.​

For each pair, it plots the score as a function of the bin centre (0.5, 1.5, …, 19.5 Å).​

Plots are saved as PNG or PDF and give a quick view of preferred and disfavoured distances.​

3. Scoring new RNA structures
The script score_structure.py applies the potential to a given RNA 3D model.​

Steps:

Parse the PDB file and extract C3' atoms for each chain.​

For each valid residue pair (intra‑chain, sequence separation ≥ 3), compute the Euclidean distance between C3' atoms.​

Identify the corresponding nucleotide pair type (e.g. AU or UA, treated consistently).​

Obtain a score for the exact distance by linear interpolation between the two nearest bins in the pre‑computed profile.​

Sum all pair scores to obtain an estimated Gibbs free energy for the conformation; lower scores are expected for more native‑like structures.​

The script can be extended to output per‑residue or per‑pair contributions to help diagnose misfolded regions.​

Command‑line usage (planned)
Examples of intended usage are:

bash
# 1) Train the potential
python scripts/training.py \
    --pdb_dir data/train_pdb \
    --max_dist 20.0 \
    --bin_width 1.0 \
    --out_dir scores/

# 2) Plot interaction profiles
python scripts/plot_profiles.py \
    --scores_dir scores/ \
    --out_dir plots/

# 3) Score a candidate RNA structure
python scripts/score_structure.py \
    --pdb data/test_pdb/model.pdb \
    --scores_dir scores/
These commands may be adapted depending on the actual script implementation and the location of the datasets.​

Implementation notes
PDB parsing and distance calculations are implemented with standard Python libraries and/or Biopython, focusing on clarity and reproducibility.​

The code is structured into small reusable functions (for example: reading PDBs, building histograms, interpolating scores) to limit redundancy and simplify testing.​

A short --help option is planned for each script to document arguments and default values.​

Version control is handled with git, and this GitHub repository will be used to track incremental improvements of the code and the analysis.​

