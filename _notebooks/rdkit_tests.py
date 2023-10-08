from rdkit import Chem

# Create a molecule object
mol = Chem.MolFromSmiles("CCO")  # Example molecule: Ethanol (C2H5OH)

# Add explicit hydrogens

mols = [mol]
print(f"Initial:                     {mols=}")

mol_with_Hs = Chem.AddHs(mol)
print(f"After Chem.AddHs(mol):       {mols=}")

mols[0] = mol_with_Hs
print(f"After mols[0] = mol_with_Hs: {mols=}")

# Iterate through the atoms in the molecule
for atom in mol.GetAtoms():
    print(f"Atom Symbol: {atom.GetSymbol()}")
