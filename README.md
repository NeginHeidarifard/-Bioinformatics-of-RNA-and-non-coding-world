
# Statistical scoring of RNA 3D structures

This project is part of the *RNA structure bioinformatics* practical course  
(Master 2 GENIOMHE, Université Paris-Saclay).

The goal is to implement a simple **statistical potential** that approximates the
Gibbs free energy of RNA 3D conformations based on **C3′–C3′ distance
distributions** between nucleotides.

> **Idea:** native RNA structures are expected to have “favourable” patterns of
> C3′–C3′ distances. We learn these patterns from solved structures and then
> use them as an **objective function** to score new models.

---

## Project goals

For a given ribonucleotide chain, the RNA folding problem consists in finding
the *native fold* among an astronomically large number of possible
conformations. The native fold is (in reality) the one with the **lowest Gibbs
free energy ΔG**. We therefore want to build an **objective function** that
estimates ΔG from a 3D structure.

In this TP, the objective function is a **pairwise statistical potential**
defined on distances between C3′ atoms:

- It is **trained** on experimentally determined RNA structures.
- It is then used to **score** new conformations, including models from the
  RNA-Puzzles dataset.

Concretely, the project:

- computes distance distributions between C3′ atoms in native RNA structures;
- derives a pairwise statistical potential for the ten nucleotide pairs  
  (AA, AU, AC, AG, UU, UC, UG, CC, CG, GG);
- uses this potential to estimate the Gibbs free energy of a candidate RNA
  conformation.

---

## Repository structure

```text
.
├── data/
│   ├── train_pdb/          # training RNA PDB files (native structures)
│   └── test_pdb/           # RNA-Puzzles or other test structures
├── scores/                 # learned scoring profiles for each base pair
├── plots/                  # interaction profiles and histograms/KDE
├── scripts/
│   ├── training.py         # Script 1: learn statistical potential
│   ├── plot_profiles.py    # Script 2: plot score profiles
│   ├── plot_hist_kde.py    # Extra: histogram + KDE for distances
│   └── score_structure.py  # Script 3: score new structures
├── README.md
└── requirements.txt

The data/ folders are not tracked by default (only local copies);
 paths can be changed via command-line arguments.

Dependencies and installation
The implementation is deliberately lightweight and uses only standard Python
 tools.
Python: ≥ 3.9 (tested with Miniconda Python 3.13 on WSL / Ubuntu)


Packages:


numpy


biopython


matplotlib


The requirements.txt file lists these dependencies:
numpy
biopython
matplotlib

Installation with conda (recommended)
# 1) Create and activate an environment
conda create -n rna_potential python=3.11 -y
conda activate rna_potential

# 2) Install dependencies
pip install -r requirements.txt


Training data
For this TP, the statistical potential is trained on a small set of known
 RNA 3D structures downloaded from the PDB:
1EHZ — hairpin ribozyme (native structure)


4TNA, 6TNA — tRNA^Asp structures


These files are stored in:
data/train_pdb/
    1EHZ.pdb
    4TNA.pdb
    6TNA.pdb

For evaluation, several RNA-Puzzles models were downloaded and placed in:
data/test_pdb/
    1Y26.pdb
    2L8H.pdb
    5T5A.pdb
    ...


Method overview
1. Training the statistical potential (training.py)
The training script implements exactly what is described in the TP:
“Compute interatomic distances from the given dataset of PDB files;
 only C3′ atoms are taken into account;
 10 distance distributions for the 10 base pairs;
 only intrachain distances;
 only residues separated by at least 3 positions (i, i+4, i+5, …);
 compute observed and reference frequencies;
 compute the log-ratio; cap maximum score at 10;
 generate 10 files of 20 lines (one file per base pair).”
Atom and pair selection
Only C3′ atoms are used.


C3′ is in the backbone, always present even when bases are partially
 unresolved.


Only intra-chain distances are considered.


Only nucleotide types A, U, C, G are kept.


Only pairs with sequence separation ≥ 3 are considered:


skip pairs (i, i+1), (i, i+2), (i, i+3)


keep (i, i+4), (i, i+5), …


This removes very local geometry (almost rigid backbone) so that the model
 focuses on distances informative about 3D folding rather than trivial
 covalent constraints.
Distance binning
Distances between C3′ atoms are computed using the Euclidean norm and
 binned into 20 intervals between 0 and 20 Å:
bin width Δr = 1.0 Å


bin centres: 0.5, 1.5, …, 19.5 Å


For each of the 10 unordered base pairs:
AA, AU, AC, AG, UU, UC, UG, CC, CG, GG

we compute a histogram of C3′–C3′ distances.
Observed frequency
For a given pair type (e.g. AU) and bin r:
( N_{ij}(r) ) = number of AU pairs whose distance falls in bin r


( N_{ij} ) = total number of AU pairs over all bins


The observed probability is:
[
 f^{\mathrm{OBS}}{ij}(r) = N{ij}(r) / N_{ij}.
 ]
Reference frequency (“XX”)
The reference distribution ignores nucleotide identity:
Count all pairs of bases (A/U/C/G) together as “XX”.


( N_{XX}(r) ) = number of any base pair whose distance falls in bin r


( N_{XX} ) = total number of such pairs


[
 f^{\mathrm{REF}}{XX}(r) = N{XX}(r) / N_{XX}.
 ]
This corresponds to the “XX pair” in the TP handout.
Pseudo-energy (score)
For each base pair (i,j) and bin r the statistical potential is
[
 u_{ij}(r)
 = -\log \left( \frac{f^{\mathrm{OBS}}{ij}(r) + \epsilon}{f^{\mathrm{REF}}{XX}(r) + \epsilon} \right)
 ]
where a small pseudocount ε is added to avoid division by zero.
If a distance is seen more often than in the reference distribution,
 the ratio > 1, log > 0, and −log < 0 → favourable (low energy).


If it is seen less often, the energy is positive (unfavourable).


As required in the TP, scores are clipped to a maximum of +10.
Output score files
The training script produces 10 text files in scores/:
scores/
  AA.txt
  AU.txt
  AC.txt
  AG.txt
  UU.txt
  UC.txt
  UG.txt
  CC.txt
  CG.txt
  GG.txt

Each file contains 20 lines, one per distance bin, with a single
 floating-point value (pseudo-energy) per line.
Command-line usage
python scripts/training.py \
    --pdb_dir data/train_pdb \
    --out_dir scores \
    --max_dist 20 \
    --bin_width 1 \
    --min_sep 3

The script prints the processed PDB files and confirms that score files were
 written to scores/.

2. Plotting interaction profiles (plot_profiles.py)
The second script corresponds to:
“The second script will plot the interaction profiles: the score as a
 function of the distance.”
For each base pair, it:
reads the corresponding scores/XX.txt file;


reconstructs the distance axis from max_dist and bin_width;


plots pseudo-energy vs distance (one curve per base pair);


saves a PNG image in the plots/ folder.


Command-line usage:
python scripts/plot_profiles.py \
    --scores_dir scores \
    --out_dir plots \
    --max_dist 20 \
    --bin_width 1

Result: 10 PNG files such as plots/AU.png, plots/CG.png, …
These plots show which distances are preferred (negative scores) or
 disfavoured (positive scores) for each nucleotide pair.

3. Histogram + KDE analysis (plot_hist_kde.py, optional)
During the course, special attention was given to histograms and
 kernel density estimation (KDE) for a specific base pair, typically AU.
The extra script plot_hist_kde.py:
recomputes all C3′–C3′ distances from data/train_pdb/;


filters them for a given base pair (e.g. AU) or for XX (all bases pooled);


builds a histogram with a given bin width;


overlays a simple Gaussian KDE with a chosen bandwidth;


saves the figure (e.g. plots/hist_kde_AU_auto.png).


This is useful to explore how the choice of bin width and bandwidth
 affects the smoothness of the profile and to justify the choice Δr = 1 Å.
Example usage:
# AU-specific histogram + KDE
python scripts/plot_hist_kde.py \
    --pdb_dir data/train_pdb \
    --out_path plots/hist_kde_AU_auto.png \
    --pair AU \
    --max_dist 20 \
    --bin_width 1 \
    --bandwidth 1.0

# Reference “XX” distribution
python scripts/plot_hist_kde.py \
    --pdb_dir data/train_pdb \
    --out_path plots/hist_kde_XX_auto.png \
    --pair XX \
    --max_dist 20 \
    --bin_width 1 \
    --bandwidth 1.0


4. Scoring new structures (score_structure.py)
The third script implements:
“Compute all the distances for a given structure (same thresholds: 20 Å
 and i, i+4); for each distance, compute a scoring value by linear
 interpolation; by summing all these scores, calculate the estimated Gibbs
 free energy of the evaluated RNA conformation.”
Given:
a PDB file (candidate structure),


the folder of pre-computed scores (scores/),


the script:
Parses the PDB file with Biopython’s PDBParser.


Extracts C3′ atoms for all valid residues (A, U, C, G).


Considers all intra-chain pairs with separation ≥ min_sep.


Computes their C3′–C3′ distances.


Identifies the base-pair type (AA, AU, …, GG).


Uses linear interpolation between the two nearest bins in the
 corresponding score profile to obtain a pseudo-energy value.


Sums the contributions of all pairs to obtain the total score.


Distances > max_dist are ignored. If a base pair type has no profile,
 the pair is skipped.
Command-line usage:
python scripts/score_structure.py \
    --pdb data/test_pdb/1Y26.pdb \
    --scores_dir scores \
    --max_dist 20 \
    --bin_width 1 \
    --min_sep 3

The script reports:
number of contributing pairs;


the parameters used (max_dist, bin_width, min_sep);


the estimated pseudo-Gibbs free energy (arbitrary units).



Full workflow
Train the potential

 python scripts/training.py \
    --pdb_dir data/train_pdb \
    --out_dir scores \
    --max_dist 20 \
    --bin_width 1 \
    --min_sep 3


Visualise interaction profiles

 python scripts/plot_profiles.py \
    --scores_dir scores \
    --out_dir plots \
    --max_dist 20 \
    --bin_width 1


(Optional) Inspect histograms + KDE

 python scripts/plot_hist_kde.py \
    --pdb_dir data/train_pdb \
    --out_path plots/hist_kde_AU_auto.png \
    --pair AU \
    --max_dist 20 \
    --bin_width 1 \
    --bandwidth 1.0


Score RNA-Puzzles models

 python scripts/score_structure.py \
    --pdb data/test_pdb/1Y26.pdb \
    --scores_dir scores \
    --max_dist 20 \
    --bin_width 1 \
    --min_sep 3

python scripts/score_structure.py \
    --pdb data/test_pdb/2L8H.pdb \
    --scores_dir scores \
    --max_dist 20 \
    --bin_width 1 \
    --min_sep 3

python scripts/score_structure.py \
    --pdb data/test_pdb/5T5A.pdb \
    --scores_dir scores \
    --max_dist 20 \
    --bin_width 1 \
    --min_sep 3



Example results on RNA-Puzzles models
Using the potential trained on 1EHZ/4TNA/6TNA, the following scores were
 obtained (energies in arbitrary units):
Structure
Role
# C3′–C3′ pairs
Total pseudo-ΔG (a.u.)
1EHZ.pdb
Native (training set)
357
~ 4.4
1Y26.pdb
RNA-Puzzles model
504
~ 135.5
2L8H.pdb
RNA-Puzzles model
155
~ 99.2
5T5A.pdb
RNA-Puzzles model
506
~ 231.5

Even though the potential is extremely simple and trained on only three
 structures, it already shows the expected trend:
the native structure 1EHZ has a much lower pseudo-energy;


RNA-Puzzles models are penalised with higher scores, sometimes by
 ~30–50×.


This demonstrates that the statistical potential discriminates between a
 known native state and less native-like predictions.

Design choices and discussion
C3′ atoms only.
 Using C3′ focuses on backbone geometry and avoids missing atoms in bases.
 It also makes the potential robust across different resolutions.


Minimum sequence separation (min_sep = 3).
 This removes trivial, local contacts determined by covalent geometry and
 keeps long-range interactions that reflect the actual fold.


Bin width (Δr = 1 Å).
 Smaller bins would be too noisy given the small dataset; larger bins would
 smear interesting peaks. Δr = 1 Å was a good compromise after inspecting
 histograms and KDEs.


Reference “XX” distribution.
 Normalising by XX makes the potential relative: we do not care
 whether a distance is frequent in absolute terms, but whether it is more
 or less frequent for this particular base pair than for generic pairs.


Pseudo-energy, not real ΔG.
 The scores are in arbitrary units, and nothing guarantees that the
 lowest score corresponds exactly to the physical Gibbs free energy.
 However, as an objective function for ranking conformations, this is
 already useful.



Limitations and possible extensions
The training dataset is very small (three PDB files); real
 statistical potentials are trained on hundreds or thousands of structures.


Only one atom per nucleotide is used; including several atoms or
 base-specific interaction centres would give richer potentials.


Entropic effects and solvent are ignored.


The potential is entirely pairwise; no multi-body or angular terms.


Possible extensions (beyond the TP):
learn separate potentials for different RNA environments (helices, loops);


add base-stacking–specific terms;


use more sophisticated KDE-based potentials instead of coarse binning.



Reproducibility and versioning
The project is tracked with git and hosted on GitHub in the repository
 RNA-Statistical-Potential.


All scripts are self-contained, documented with --help options, and
 rely only on the packages listed in requirements.txt.


The workflow described above (training → plotting → scoring) can be
 reproduced on any machine with Python ≥ 3.9 by following the exact
 command-line calls.


Author: Negin Heidarifard
 Course: RNA structure bioinformatics – Master 2 GENIOMHE, Université Paris-Saclay


