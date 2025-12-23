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
