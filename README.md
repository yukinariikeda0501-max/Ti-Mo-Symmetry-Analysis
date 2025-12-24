# Ti-Mo Phase Diagram Analysis Tools

This repository contains the supplemental materials for the paper:

**"Calculation of the Ti–Mo phase diagram using density functional theory and crystal symmetry"** by Y. Ikeda and A. Ishii.

**Overview** These codes are developed for the calculation of cluster configurations and symmetry properties in Ti-Mo alloys using Density Functional Theory (DFT) and the Cluster Variation Method (CVM). The programs detect user-specified clusters within an input atomic model and output detailed cluster information, including the specific atomic arrangement and the weighting based on the number of unit cells sharing the cluster.

**Usage** A **POSCAR**-formatted atomic model is required as input for the configuration counters. Additionally, the lattice constant and the indices of the atoms constituting the cluster must be specified via command-line arguments.

**Disclaimer** Please note that these source codes were originally developed for personal research use and may not be fully optimized for user-friendliness. Please feel free to contact the authors if you have any questions.

---

## Repository Contents

This repository consists of three main sets of tools:
1.  **Symmetry Analysis Suite**: Calculates invariant symmetry operations for atomic clusters.
2.  **Cluster Configuration Counter**: Calculates the frequency of specific atomic cluster configurations ($n_l$) in supercells.
3.  **Configuration Multiplicity Calculator**: Calculates the multiplicity ($\alpha_l$) of atomic arrangements.

All tools use a consistent floating-point tolerance (`EPS 1e-2`) and require the standard math library (`-lm`) for compilation.

---

## 1. Symmetry Analysis Suite

These tools calculate the number of invariant symmetry operations for user-defined atomic clusters. The outputs correspond to variables in **Equation (15)** and **Equation (17)** of the referenced paper.

### Theoretical Background

#### Symmetry of Cluster A ($|\mathcal{N}_{A}(G)|$)
The program outputs "Total valid operations," corresponding to $|\mathcal{N}_{A}(G)|$. This value is used to calculate the multiplicity $\mu_A$ of the cluster (Eq. 15):

$$
\mu_{A} = \frac{|\Gamma_{\phi}|}{|\mathcal{N}_{A}(G)|}
$$

(Where $|\Gamma_{\phi}|$ is `NUM_OPE`, i.e., 24 or 48 depending on the phase).

#### Intersection Group
The "Re-check" function counts the number of operations that remain valid when the rotation center is shifted to a subcluster $A'$. This corresponds to the order of the intersection group $|\mathcal{N}_{A}(G) \cap \mathcal{N}_{A'}(G)|$. This is used to calculate the number of subclusters $M(A', A)$ (Eq. 17):

$$
M(A', A) = \sum \frac{|\mathcal{N}_{A}(G)|}{|\mathcal{N}_{A}(G) \cap \mathcal{N}_{A'}(G)|}
$$

### Included Tools & Usage

Compile using `gcc` linking the math library (`-lm`).

* **BCC Phase** ($O_h$, 48 Ops, 9 Ref Atoms)
    * Source: `M_beta_output.c`
    * Compile: `gcc -o calc_beta Symmetry-Analysis/M_beta_output.c -lm`

* **HCP Phase** ($D_{6h}$, 24 Ops, 18 Ref Atoms)
    * Source: `M_alpha_output.c`
    * Compile: `gcc -o calc_alpha Symmetry-Analysis/M_alpha_output.c -lm`

* **Omega Phase** (24 Ops, 22 Ref Atoms)
    * Source: `M_omega_output.c`
    * Compile: `gcc -o calc_omega Symmetry-Analysis/M_omega_output.c -lm`

---

## 2. Cluster Configuration Counter ($n_l$ Calculator)

These programs calculate **$n_l(A)$**, which is the number of clusters appearing in the atomic models (supercells) that correspond to a specific equivalent arrangement $l$.

### Theoretical Background

This value is directly used in **Equation (19)** of the paper to calculate the probability $x_l(A)$, which is essential for determining the entropy factor $s(A)$ in the Cluster Variation Method (CVM).

$$
x_{l}(A) = \frac{\mu_{\phi}}{N\mu_{A}}n_{l}(A)
$$

(Where $N$ is the number of atoms in the model, and $\mu$ represents multiplicity).

### Included Tools

The programs are categorized by the crystal phase and the size of the clusters.

#### 1. BCC Phase ($O_h$)
* `cluster_identify_2atom_bcc.c`: Calculates $n_l$ for 2-atom clusters.
* `cluster_identify_3atom_bcc.c`: Calculates $n_l$ for 3-atom clusters.
* `cluster_identify_4atom_bcc.c`: Calculates $n_l$ for 4-atom clusters.
* `cluster_identify_5atom_bcc.c`: Calculates $n_l$ for 5-atom clusters.

#### 2. Omega Phase ($D_{6h}$)
* `cluster_identify_2atom_omega.c`: Calculates $n_l$ for 2-atom clusters.
* `cluster_identify_3atom_omega.c`: Calculates $n_l$ for 3-atom clusters.
* `cluster_identify_4atom_omega.c`: Calculates $n_l$ for 4-atom clusters.
* `cluster_identify_5atom_omega.c`: Calculates $n_l$ for 5-atom clusters.

#### 3. HCP Phase ($D_{6h}$)
* `cluster_identify_2atom_hcp.c`: Calculates $n_l$ for 2-atom clusters.
* `cluster_identify_3atom_hcp.c`: Calculates $n_l$ for 3-atom clusters.
* `cluster_identify_4atom_hcp.c`: Calculates $n_l$ for 4-atom clusters.
* `cluster_identify_5atom_hcp.c`: Calculates $n_l$ for 5-atom clusters.

### Usage

Use `gcc` (or any standard C compiler) and link the math library (`-lm`).

```bash
# Example: Compile the 3-atom counter for Omega phase
gcc -o count_omega_3 cluster_identify/omega/cluster_identify_3atom_omega.c -lm
