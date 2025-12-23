# Ti-Mo Phase Diagram Symmetry Analysis Suite

This repository contains a suite of C programs designed to analyze the crystal symmetry of Titanium-Molybdenum (Ti-Mo) alloys. These tools calculate the number of invariant symmetry operations for user-defined atomic clusters.

This software was developed for the research paper:

**"Calculation of the Ti–Mo phase diagram using density functional theory and crystal symmetry"**

## Theoretical Background

The outputs of these programs directly correspond to the variables required for **Equation (15)** and **Equation (17)** in the referenced paper.

### 1. Symmetry of Cluster A ($|\mathcal{N}_{A}(G)|$):

The program outputs "Total valid operations," which corresponds to 
$|\mathcal{N}_{A}(G)|$
. 
This value is used to calculate the multiplicity 
$\mu_A$
of the cluster (Eq. 15):

$$
\mu_{A} = \frac{|\Gamma_{\phi}|}{|\mathcal{N}_{A}(G)|}
$$

(Where 
$|\Gamma_{\phi}|$
is `NUM_OPE`, i.e., 24 or 48 depending on the phase).

### 2. Intersection Group

The "Re-check" function counts the number of operations that remain valid when the rotation center is shifted to a subcluster $A'$. This corresponds to the order of the intersection group:

$$
|\mathcal{N}_{A}(G) \cap \mathcal{N}_{A'}(G)|
$$

This is used to calculate the number of subclusters $M(A', A)$ (Eq. 17):

$$
M(A', A) = \sum \frac{|\mathcal{N}_{A}(G)|}{|\mathcal{N}_{A}(G) \cap \mathcal{N}_{A'}(G)|}
$$

## Included Tools

All tools use a consistent floating-point tolerance (`EPS 1e-2`).

* **M_bcc_output.c (BCC Phase):** 48 Symmetry Operations ($O_h$), 9 Reference Atoms.
* **M_hcp_output.c (HCP Phase):** 24 Symmetry Operations ($D_{6h}$), 18 Reference Atoms.
* **M_w_output.c (Omega Phase):** 24 Symmetry Operations, 22 Reference Atoms.

## Usage

### Compilation

Use `gcc` (or any standard C compiler) and link the math library (`-lm`).

```bash
gcc -o calc_bcc M_bcc_output.c -lm
gcc -o calc_hcp M_hcp_output.c -lm
gcc -o calc_w   M_w_output.c   -lm


# Ti-Mo Cluster Configuration Counter ($n_l$ Calculator)

This repository contains a suite of C programs designed to calculate the frequency of specific atomic cluster configurations within Titanium-Molybdenum (Ti-Mo) supercells.

This software was developed for the research paper:

**"Calculation of the Ti–Mo phase diagram using density functional theory and crystal symmetry"**

## Overview

These programs calculate **$n_l(A)$**, which is the number of clusters appearing in the atomic models (supercells) that correspond to a specific equivalent arrangement $l$.

This value is directly used in **Equation (19)** of the paper to calculate the probability $x_l(A)$, which is essential for determining the entropy factor $s(A)$ in the Cluster Variation Method (CVM).

$$
x_{l}(A) = \frac{\mu_{\phi}}{N\mu_{A}}n_{l}(A)
$$

(Where $N$ is the number of atoms in the model, and $\mu$ represents multiplicity).

## Included Tools

The programs are categorized by the crystal phase and the size of the clusters.
*(Note: Tools for HCP 4-atom and 5-atom clusters also exist and follow the same logic, though they may not be listed below.)*

### 1. BCC Phase ($O_h$)
* `cluster_idenfity_2atom_bcc.c`: Calculates $n_l$ for 2-atom clusters.
* `cluster_idenfity_3atom_bcc.c`: Calculates $n_l$ for 3-atom clusters.
* `cluster_idenfity_4atom_bcc.c`: Calculates $n_l$ for 4-atom clusters.
* `cluster_idenfity_5atom_bcc.c`: Calculates $n_l$ for 5-atom clusters.

### 2. Omega Phase ($D_{6h}$)
* `cluster_idenfity_2atom_omega.c`: Calculates $n_l$ for 2-atom clusters.
* `cluster_idenfity_3atom_omega.c`: Calculates $n_l$ for 3-atom clusters.
* `cluster_idenfity_4atom_omega.c`: Calculates $n_l$ for 4-atom clusters.
* `cluster_idenfity_5atom_omega.c`: Calculates $n_l$ for 5-atom clusters.

### 3. HCP Phase ($D_{6h}$)
* `cluster_idenfity_2atom_hcp.c`: Calculates $n_l$ for 2-atom clusters.
* `cluster_idenfity_3atom_hcp.c`: Calculates $n_l$ for 3-atom clusters.
* *(4-atom and 5-atom counters are also part of this suite)*

## Usage

### Compilation

Use `gcc` (or any standard C compiler) and link the math library (`-lm`).

```bash
# Example: Compile the 3-atom counter for Omega phase
gcc -o count_omega_3 cluster_idenfity_3atom_omega.c -lm
