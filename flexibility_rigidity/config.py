# Atom-specific radii
ATOM_RADII = {
    "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
    "F": 1.47, "Cl": 1.75, "Br": 1.85, "I": 1.98
}

# Supported ligand elements
LIGAND_ELEMENTS = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]

# Feature index mapping for atom pairs
FEATURE_INDEX = {
    "CC": 0, "CN": 1, "CO": 2, "CS": 3, "CP": 4, "CF": 5, "CCl": 6, "CB": 7, "CI": 8,
    "NC": 9, "NN": 10, "NO": 11, "NS": 12, "NP": 13, "NF": 14, "NCl": 15, "NB": 16, "NI": 17,
    "OC": 18, "ON": 19, "OO": 20, "OS": 21, "OP": 22, "OF": 23, "OCl": 24, "OB": 25, "OI": 26,
    "SC": 27, "SN": 28, "SO": 29, "SS": 30, "SP": 31, "SF": 32, "SCl": 33, "SB": 34, "SI": 35
}

# Hyperparameters for kernel computation
HYPERPARAMS = {
    "eta_factors": [1, 12.5],
    "lorentzian_powers": [40, 40],
    "distance_thresholds": [39, 47]
}
