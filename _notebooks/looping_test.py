import numpy as np
from rdkit import Chem


def assign_isotopes(arr, prefix):
    if isinstance(arr, Chem.Mol):
        print("".join([f"[{i}]" for i in prefix]) + Chem.MolToSmiles(arr))
    else:
        for i in range(len(arr)):
            new_prefix = prefix + [i]
            arr[i] = assign_isotopes(arr[i], new_prefix)
    return arr

# Example usage with a 3-dimensional ndarray containing RDKit molecule objects
mol1 = Chem.MolFromSmiles("CCO")
mol2 = Chem.MolFromSmiles("CCN")
mol3 = Chem.MolFromSmiles("CCOCC")
mol4 = Chem.MolFromSmiles("CCCN")
arr = np.array([[[mol1, mol2], [mol3, mol4]], [[mol2, mol3], [mol4, mol1]]], dtype=object)

prefix = []  # Initialize the prefix as an empty list
arr = assign_isotopes(arr, prefix)
print("After looping:")
print(arr)