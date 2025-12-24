# Ti-Mo Phase Diagram Analysis Tools

This repository contains the supplemental materials for the paper:

**"Calculation of the Ti–Mo phase diagram using density functional theory and crystal symmetry"** by Y. Ikeda and A. Ishii.

**Overview**
This repository provides a comprehensive suite of tools for analyzing Ti-Mo alloys using Density Functional Theory (DFT) and the Cluster Variation Method (CVM). The tools primarily calculate:
1.  **Symmetry Properties:** Analyzes invariant symmetry operations for clusters.
2.  **$n_l$ (Cluster Configurations):** Detects user-specified clusters within an input atomic model and calculates their frequency, considering the number of unit cells sharing the cluster.
3.  **$\alpha_l$ (Multiplicity):** Calculates the number of equivalent atomic arrangements for specific clusters.

**Usage**
Usage depends on the specific tool set:
* **Cluster Configuration Counters ($n_l$):** Require a **POSCAR**-formatted atomic model as input. The lattice constant and atom indices must be specified via command-line arguments.
* **Symmetry & Multiplicity Tools:** Generally operate with interactive inputs or pre-defined parameters.
* *Please refer to the specific sections below for detailed instructions on each tool.*

**Disclaimer**
Please note that these source codes were originally developed for personal research use and may not be fully optimized for user-friendliness. Please feel free to contact the authors if you have any questions.

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
```

---

## 3. Configuration Multiplicity Calculator ($\alpha_l$ Calculator)

These programs calculate the multiplicity **$\alpha_l(A)$**, which represents the number of equivalent atomic arrangements for a specific configuration $l$ on a given cluster $A$.

### Theoretical Background

This value appears in **Equation (18)** of the referenced paper and is crucial for calculating the entropy factor $s(A)$.

$$
s(A) = k_B \sum_{l} \alpha_{l}(A) x_{l}(A) \ln x_{l}(A)
$$

* **$\alpha_l(A)$**: The number of arrangements equivalent to arrangement $l$ under the symmetry operations of the cluster, $\mathcal{N}_A(G)$.
* **Orbits**: The program classifies all $2^N$ possible binary configurations (Ti/Mo) into "orbits" (groups of equivalent configurations). The size of each orbit corresponds to $\alpha_l$.
* **Verification**: This tool automates the calculation of values such as those shown in **Figure 7** for the $\beta_5$ cluster.

### Included Tools

These tools take a specific set of atom indices (defining a cluster) as input and output the multiplicity for all unique configurations.

* **BCC Phase** ($O_h$ Symmetry)
    * Source: `beta_a_l.c`
    * Compile: `gcc -o beta_a_l beta_a_l.c -lm`
* **HCP Phase** ($D_{6h}$ Symmetry)
    * Source: `alpha_a_l.c`
    * Compile: `gcc -o alpha_a_l alpha_a_l.c -lm`
* **Omega Phase** ($D_{6h}$ Symmetry)
    * Source: `omega_a_l.c`
    * Compile: `gcc -o omega_a_l omega_a_l.c -lm`

### Usage

1.  **Compile** the program using `gcc` with the math library linked (`-lm`).
2.  **Run** the executable.
3.  **Input** the number of atoms in the cluster.
4.  **Input** the specific atom indices defining the cluster (based on the reference atomic positions in the code).

#### Example: Calculating $\alpha_l$ for $\beta_5$ cluster (Figure 7 in paper) using `beta_a_l`

**Input:**
```text
Atoms: 3
Atom Indices: 1, 2, 5
```

**Output (Conceptual):**
The program will generate a table classifying configurations by grouping equivalent Ti/Mo arrangements.

| ID | Ti Atoms | Mo Atoms | Alpha_l |
| :--- | :--- | :--- | :--- |
| 1 | 1,2,5 | (None) | 1 |
| 2 | 5 | 1,2 | 2 |
| ... | ... | ... | ... |

* **ID**: Unique identifier for the configuration orbit.
* **Ti/Mo Atoms**: Shows which atoms are occupied by Ti or Mo for a representative configuration of that orbit.
* **Alpha_l**: The multiplicity (degeneracy) of that configuration.
