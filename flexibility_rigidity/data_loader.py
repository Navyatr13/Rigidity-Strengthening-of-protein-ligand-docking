import os
import numpy as np
import pandas as pd
from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb

def load_datasets():
    """
    Load training and test datasets.
    """
    train_data = pd.read_csv("v2007_refine_listt.csv")
    test_data = pd.read_csv("v2007_core_listt.csv")
    train_list = np.asarray(train_data)[:, 0]
    test_list = np.asarray(test_data)[:, 0]
    return train_list, test_list


def load_molecule_and_protein(file_path, ligand_file, protein_file):
    """
    Load a ligand and protein from the specified files.
    """
    pmol = PandasMol2().read_mol2(os.path.join(file_path, ligand_file))
    mol_list = pmol.df[pmol.df["atom_type"] != "H"][["atom_name", "x", "y", "z"]]

    ppdb = PandasPdb().read_pdb(os.path.join(file_path, protein_file))
    prot_list = ppdb.df["ATOM"][ppdb.df["ATOM"]["atom_name"] != "H"][["atom_name", "x_coord", "y_coord", "z_coord"]]

    return np.asarray(mol_list), np.asarray(prot_list)
