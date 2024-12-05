# Flexibility Rigidity Index with Multiple Kernels

This script calculates the **Flexibility Rigidity Index** using multiple kernels and different hyperparameters for ligand-protein interactions. The code processes `.pdb` and `.mol2` files, computes interaction distances, and applies exponential and Lorentzian kernels based on atom types and their distances.

---

## Features

- **Input Files**:
  - Reads `.pdb` and `.mol2` files for protein and ligand information.
  - Supports dataset lists in `.csv` format for train and test sets.
- **Atom-Specific Kernels**:
  - Computes Euclidean distances between atoms.
  - Applies exponential and Lorentzian kernels based on ligand and protein atom types.
- **Output**:
  - Generates a feature matrix saved as a CSV file.
  - Feature matrix includes flexibility rigidity scores for specified atom pairs.

---

## Prerequisites

### Libraries
- `numpy`: For numerical computations.
- `pandas`: For handling `.csv` data.
- `biopandas`: For reading `.pdb` and `.mol2` files.
- `os`, `itertools`, `math`, `sys`: For file handling and mathematical operations.

### Installation
Install required packages using pip:
``` pip install numpy pandas biopandas ```


---

## File Structure
- **Input Data**:
  - `/v2007/`: Directory containing `.pdb` and `.mol2` files.
  - `v2007_refine_listt.csv`: CSV containing the training data file names.
  - `v2007_core_listt.csv`: CSV containing the test data file names.
- **Script**:
  - `flexibility_rigidity.py`: Main script for computing the flexibility rigidity index.
- **Output**:
  - `Lor_2kerneltrain_2.5-1-12_2.5-5.5-21.csv`: CSV file containing the computed feature matrix.

---

## Usage

### Step 1: Prepare Input Data
1. Place `.pdb` and `.mol2` files in the `/v2007/` directory.
2. Create train and test CSV files:
   - `v2007_refine_listt.csv` and `v2007_core_listt.csv`.
   - Each file should contain file names (one per line) without extensions.

### Step 2: Run the Script
Execute the script:
``` python flexibility_rigidity.py ```


### Step 3: Output
The script generates a CSV file (`Lor_2kerneltrain_2.5-1-12_2.5-5.5-21.csv`) with the feature matrix.

---

## Code Overview

### Key Functions
1. **Atom Distance Computation**:
   - Calculates Euclidean distance between protein and ligand atoms.
   - Filters distances below a specified threshold.
2. **Kernel Application**:
   - Applies exponential and Lorentzian kernels using atom-specific radii.
   - Uses hyperparameters (e.g., `t1`, `k1`, `v1`) to compute kernel values.
3. **Feature Matrix Construction**:
   - Builds a matrix with flexibility rigidity scores for specific atom pairs.

---

## Hyperparameters

| Parameter | Description                               | Default Value |
|-----------|-------------------------------------------|---------------|
| `t1`      | Scaling factor for exponential kernel 1   | `1`           |
| `t2`      | Scaling factor for exponential kernel 2   | `12.5`        |
| `k1`      | Power for kernel 1                        | `5`           |
| `k2`      | Power for kernel 2                        | `5`           |
| `v1`      | Parameter for Lorentzian kernel 1         | `40`          |
| `v2`      | Parameter for Lorentzian kernel 2         | `40`          |

---

## Example Dataset

Example training data in `v2007_refine_listt.csv`:
``` file1 file2 file3 ... ```


Each file should correspond to `.pdb` and `.mol2` files in the `/v2007/` directory:
- `file1_protein.pdb`
- `file1_ligand.mol2`

---

## Future Enhancements

1. Extend kernel calculations to include additional molecular properties.
2. Add visualization for kernel values and distance distributions.
3. Automate dataset preparation and validation.

---

