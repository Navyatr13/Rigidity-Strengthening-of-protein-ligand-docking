import numpy as np
import itertools
from config import LIGAND_ELEMENTS, ATOM_RADII, FEATURE_INDEX, HYPERPARAMS


def calculate_features(train_list):
    """
    Calculate feature matrix for the provided training list.
    """
    num_samples = len(train_list)
    num_features = len(FEATURE_INDEX) * 2  # 36 features for two kernels
    feature_matrix = np.zeros((num_samples, num_features))

    for i, file_name in enumerate(train_list):
        mol_list, prot_list = load_molecule_and_protein(
            "/mnt/home/tippanan/prj/v2007",
            f"{file_name}_ligand.mol2",
            f"{file_name}_protein.pdb"
        )
        for x, y in itertools.product(prot_list, mol_list):
            compute_feature_for_pair(feature_matrix, i, x, y, kernel_idx=0)
            compute_feature_for_pair(feature_matrix, i, x, y, kernel_idx=1)

    return feature_matrix


def compute_feature_for_pair(feature_matrix, row_idx, x, y, kernel_idx):
    """
    Compute kernel-based feature for a pair of atoms.
    """
    dist = np.linalg.norm(x[1:4] - y[1:4], axis=0)
    if dist < HYPERPARAMS["distance_thresholds"][kernel_idx]:
        r_i, r_j = ATOM_RADII[x[0][0]], ATOM_RADII[y[0][0]]
        eta = HYPERPARAMS["eta_factors"][kernel_idx] * (r_i + r_j)
        lor_v = 1 / (1 + (dist / eta) ** HYPERPARAMS["lorentzian_powers"][kernel_idx])
        feature_idx = FEATURE_INDEX[x[0][0] + y[0][0]] + kernel_idx * 36
        feature_matrix[row_idx, feature_idx] += lor_v
