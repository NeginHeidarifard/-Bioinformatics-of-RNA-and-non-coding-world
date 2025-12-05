
---

**Statistical scoring of RNA 3D structures**
This project is part of the *RNA structure bioinformatics* practical course
(Master 2 GENIOMHE, Université Paris-Saclay).

The goal is to implement a simple **statistical potential** that approximates the Gibbs free energy of RNA 3D conformations based on **C3′–C3′ distance distributions** between nucleotides.

> Native RNA structures exhibit characteristic, non-random C3′–C3′ distance patterns.
> We learn these patterns from solved RNA structures and use them as an objective function to score new models.

---

## **Project goals**

We aim to build an **objective function** that estimates ΔG from an RNA 3D structure.

The objective function is a **pairwise statistical potential** based on distances between C3′ atoms:

* It is **trained** on experimentally determined native RNA structures.
* It is then **used to score** new conformations, including models from the RNA-Puzzles dataset.

Concretely, the project:

* computes distance distributions between C3′ atoms in native structures,
* derives potentials for **ten nucleotide pairs**:
  *AA, AU, AC, AG, UU, UC, UG, CC, CG, GG*,
* uses these potentials to estimate a pseudo–Gibbs free energy for new candidate models.

---

## **Repository structure**

```text
.
├── data/
│   ├── train_pdb/        # native training RNA structures (1EHZ, 4TNA, 6TNA)
│   └── test_pdb/         # RNA-Puzzles or evaluation structures
├── scores/               # learned statistical potentials for each pair
├── plots/                # interaction profiles and histogram/KDE plots
├── scripts/
│   ├── training.py       # Script 1: learn statistical potential
│   ├── plot_profiles.py  # Script 2: plot score profiles
│   ├── plot_hist_kde.py  # Optional: histogram + KDE exploration
│   └── score_structure.py# Script 3: score new RNA structures
├── requirements.txt
└── README.md
```

---

## **Dependencies & installation**

This project uses only standard scientific Python tools.

### **Python version**

* Python **≥ 3.9**
* Tested on Conda Python **3.13** under WSL/Ubuntu

### **Packages**

```
numpy
biopython
matplotlib
```

### **Install via conda**

```bash
conda create -n rna_potential python=3.11 -y
conda activate rna_potential
pip install -r requirements.txt
```

---

## **Training data**

Native structures used for training:

| PDB ID | Description               |
| ------ | ------------------------- |
| 1EHZ   | Hairpin ribozyme (native) |
| 4TNA   | tRNA^Asp                  |
| 6TNA   | tRNA^Asp                  |

Stored in:

```
data/train_pdb/
```

RNA-Puzzles or other models for evaluation:

```
data/test_pdb/
    1Y26.pdb
    2L8H.pdb
    5T5A.pdb
    ...
```

---

# **Method overview**

---

## **1. Learning the statistical potential — `training.py`**

The script implements the TP specification exactly:

> * compute all C3′–C3′ distances
> * restrict to intrachain pairs
> * require sequence separation ≥ 3
> * build distance histograms for 10 base pairs
> * compute observed vs reference frequencies
> * compute −log(f_obs / f_ref)
> * cap scores at +10
> * output 20-bin profiles (0–20 Å, bin size 1 Å)

### **Atom and pair filtering**

* Only **C3′ atoms**
* Only **A, U, C, G** nucleotides
* Only **intrachain** distances
* Only pairs with **sequence separation ≥ 3**

Removes trivial local geometry (i,i+1,i+2,i+3).

### **Distance binning**

* Range: **0–20 Å**
* Bin width: **1.0 Å**
* 20 bins, centers at 0.5, 1.5, …, 19.5 Å

### **Observed vs reference frequencies**

* For each base pair type ( ij ): compute histogram ( f^{OBS}_{ij}(r) )
* Reference distribution “XX” pools all nucleotides
* Statistical potential:

[
u_{ij}(r) = -\log\left(\frac{f^{OBS}*{ij}(r) + \epsilon}{f^{REF}*{XX}(r) + \epsilon}\right)
]

with pseudocount ε.

### **Output**

10 files in `scores/`:

```
AA.txt AU.txt AC.txt AG.txt
UU.txt UC.txt UG.txt
CC.txt CG.txt GG.txt
```

Each contains **20 pseudo-energy values**.

### **Run**

```bash
python scripts/training.py \
  --pdb_dir data/train_pdb \
  --out_dir scores \
  --max_dist 20 \
  --bin_width 1 \
  --min_sep 3
```

---

## **2. Plotting interaction profiles — `plot_profiles.py`**

Generates the “score vs distance” curves for each base pair.

### Produces:

```
plots/AA.png ... plots/GG.png
```

### Run:

```bash
python scripts/plot_profiles.py \
  --scores_dir scores \
  --out_dir plots \
  --max_dist 20 \
  --bin_width 1
```

---

## **3. Histogram + KDE exploration — `plot_hist_kde.py` (optional)**

Used to justify bin width and smoothing.

* Reads C3′–C3′ distances (AU, or XX)
* Builds histogram with chosen bin width
* Computes simple Gaussian KDE
* Saves combined figure

### Example:

```bash
python scripts/plot_hist_kde.py \
  --pdb_dir data/train_pdb \
  --out_path plots/hist_kde_AU_auto.png \
  --pair AU \
  --max_dist 20 \
  --bin_width 1 \
  --bandwidth 1.0
```

---

# **4. Scoring new structures — `score_structure.py`**

Implements:

> compute all distances for a structure → interpolate → sum → estimated ΔG

### Steps:

* parse PDB
* extract C3′ atoms
* compute valid distances
* identify base-pair type
* **linear interpolation** between adjacent bins
* sum energies to obtain pseudo–Gibbs free energy

### Run example:

```bash
python scripts/score_structure.py \
  --pdb data/test_pdb/1Y26.pdb \
  --scores_dir scores \
  --max_dist 20 \
  --bin_width 1 \
  --min_sep 3
```

---

# **Full workflow**

```bash
# 1. Train
python scripts/training.py --pdb_dir data/train_pdb --out_dir scores --max_dist 20 --bin_width 1 --min_sep 3

# 2. Plot pair profiles
python scripts/plot_profiles.py --scores_dir scores --out_dir plots --max_dist 20 --bin_width 1

# 3. (Optional) histogram + KDE
python scripts/plot_hist_kde.py --pdb_dir data/train_pdb --out_path plots/hist_kde_AU_auto.png --pair AU --bin_width 1 --bandwidth 1.0

# 4. Score RNA-Puzzles models
python scripts/score_structure.py --pdb data/test_pdb/1Y26.pdb --scores_dir scores --max_dist 20 --bin_width 1 --min_sep 3
```

---

# **Example results**

| Structure | Role   | # pairs | Pseudo-ΔG |
| --------- | ------ | ------- | --------- |
| 1EHZ      | Native | 357     | ~4.4      |
| 1Y26      | Puzzle | 504     | ~135.5    |
| 2L8H      | Puzzle | 155     | ~99.2     |
| 5T5A      | Puzzle | 506     | ~231.5    |

Trend:

* native structure = **much lower** pseudo-energy
* Puzzle models = **higher** (30–50×)

---

# **Design choices**

* C3′ atoms: robust, always present, avoids missing data
* min separation = 3: removes rigid backbone geometry
* bin width = 1 Å: balance noise vs resolution
* reference XX distribution: gives relative, not absolute, likelihoods
* arbitrary units: not real ΔG but effective for ranking conformations

---

# **Limitations & extensions**

### Current limitations:

* tiny training set (3 structures)
* only backbone geometry (C3′)
* ignores entropy and solvent
* only pairwise terms, no angles or higher-order interactions

### Possible extensions:

* environment-specific potentials (helices vs loops)
* stacking potentials
* KDE-based continuous potentials

---

# **Reproducibility**

* version-controlled with Git
* scripts are self-contained and have `--help` documentation
* experiment workflow is fully described
* reproducible on any machine with Python ≥ 3.9

---

**Author:** Negin Heidarifard
**Course:** RNA Structure Bioinformatics – Master 2 GENIOMHE, Université Paris-Saclay

---

