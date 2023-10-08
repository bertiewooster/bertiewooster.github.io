# %%
# !pip install rdkit
# !pip install mendeleev

# %%
import copy
import math
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mendeleev import element, isotope
from rdkit import Chem
from rdkit.Chem import Descriptors

# %%
mol = Chem.MolFromSmiles(sml)

# %%
mol

# %%
class MolecularIsotope():
  """Store a molecule's isotopic properties."""
  def __init__(
      self,
      sml: str,
      mol: Chem.rdchem.Mol = None,
      abundance: float = None):
    self.sml = sml
    self.abundance = abundance
    if mol is not None:
      self.mol = mol
    else:
      try:
        self.mol = Chem.MolFromSmiles(sml)
      except:
        self.mol = Chem.MolFromSmiles(sml, sanitize=False)
  def update(self):
    self.mass = Descriptors.ExactMolWt(self.mol)
    self.canon_sml = Chem.MolToSmiles(mol)

# %%
def set_isotopes(arr:np.ndarray, sml: str, current_index=()):
    if isinstance(arr, np.ndarray):
        for i in range(len(arr)):
            arr[i] = set_isotopes(arr=arr[i], current_index=current_index + (i,), sml=sml)
        return arr
    else:
        # For an individual item in the ndarray,
        #   create molecule and set the isotopes of each of its atoms
        item = MolecularIsotope(sml=sml)
        mol = item.mol
        abundance = 1
        for atom_index, atom in enumerate(mol.GetAtoms()):
          atom_symbol = atom.GetSymbol()
          isotope_data_to_use = isotopes[atom_symbol][current_index[atom_index]]
          isotope_to_use = isotope_data_to_use[0]
          atom.SetIsotope(isotope_to_use)
          abundance *= isotope_data_to_use[1]
        item.update()
        item.abundance = abundance
        return item

# %%
num_atoms = len(mol.GetAtoms())

# Build iterable of number of isotopes by atom index
isotopes_by_atom_index = ()
for atom_index, atom in enumerate(mol.GetAtoms()):
  atom_symbol = atom.GetSymbol()
  isotopes_by_atom_index += (len(isotopes[atom_symbol]),)

# %%
# Create the empty numpy ndarray with the specified shape
mol_isotopes = np.empty(isotopes_by_atom_index, dtype=object)

# Set the isotopes for each atom in each molecule
mol_isotopes = set_isotopes(arr = mol_isotopes, sml = sml, current_index=())

# %%
mol_isotopes_flat = flatten_ndarray(mol_isotopes)
mol_isotopes_flat.sort(key=lambda x:x.mass)
len(mol_isotopes_flat)

# %%
mols_flat = [mol_isotope.mol for mol_isotope in mol_isotopes_flat]
abundance_flat = [mol_isotope.abundance for mol_isotope in mol_isotopes_flat]
mass_flat = [mol_isotope.mass for mol_isotope in mol_isotopes_flat]
legends = [f"{mass:.3f}: {abundance:.3f}" for mass, abundance in zip(mass_flat, abundance_flat)]
abundance_sum = sum(abundance_flat)
# Chem.Draw.MolsToGridImage(mols = mols_flat, legends = legends, molsPerRow=3)

# %% [markdown]
# It might be tempting to match molecules by mass, but a floating-point value can have some error or uncertainty. In some cases, the mass differences between different species can be quite small (for example less than 0.002 amu if the sum of the mass numbers is the same), so adding a tolerance risks lumping different molecules together. So it's best to use some discrete identifier. In this case, we use the SMILES, which contains isotope labels.

# %%
# Merge identical molecules
merged_molecular_isotpes = [mol_isotopes_flat[0]]
for molecular_isotope in mol_isotopes_flat[1:]:
  last_merged = merged_molecular_isotpes[-1]
  if Chem.MolToSmiles(molecular_isotope.mol) == Chem.MolToSmiles(last_merged.mol):
    last_merged.abundance += molecular_isotope.abundance
  else:
    merged_molecular_isotpes.append(molecular_isotope)

# %%
mols_merged_flat = [mol_isotope.mol for mol_isotope in merged_molecular_isotpes]
abundance_merged_flat = [mol_isotope.abundance for mol_isotope in merged_molecular_isotpes]
mass_merged_flat = [mol_isotope.mass for mol_isotope in merged_molecular_isotpes]
legends = [f"{mass:.9f}: {abundance:.3f}" for mass, abundance in zip(mass_merged_flat, abundance_merged_flat)]

# %%
for legend in legends:
  print(legend)

# %%
# Print masses and differences between adjacent molecules
for index, mass in enumerate(mass_merged_flat):
  if all([index > 0, index < len(mass_merged_flat) - 1]):
    print(f"{index:<4} {mass:<25} {mass_merged_flat[index-1]:<25} {mass - mass_merged_flat[index-1]}")

# %%
Chem.Draw.MolsToGridImage(mols = mols_merged_flat, legends = legends)

# %%
abundance_sum = sum(abundance_merged_flat)
print(abundance_sum)

# %%
df = pd.DataFrame({
    'Mass': mass_merged_flat,
    'Abundance': abundance_merged_flat
})

# %%
df = df.loc[df['Abundance'] >= 0.01]

# %%
# Create the scatter plot using Seaborn
sns.scatterplot(x='Mass', y='Abundance', data=df, marker='')

plt.xlabel('Mass')
plt.ylabel('Abundance')
plt.title(f'{sml} molecular isotopic distribution')

# Adjust x-axis limits to allow more space on the left and right for labels
xrange = max(df['Mass']) - min(df['Mass'])
xpad = xrange * 0.15
plt.xlim(min(df['Mass']) - xpad, max(df['Mass']) + xpad)

# Adjust y-axis limits to start at zero and allow more space on the top for labels
yrange = max(df['Abundance']) - min(df['Abundance'])
ypad = yrange * 0.1
plt.ylim(0, max(df['Abundance']) + ypad)

# Add vertical lines from each data point to the x-axis
for x, y in zip(df['Mass'], df['Abundance']):
    plt.vlines(x=x, ymin=0, ymax=y, colors='black')
    # Label the data points by mass
    plt.text(x=x, y=y, s=f'{x:.5f}', ha='center', va='bottom')

plt.show()

# %%
sml = "C=C=O"

# sml = "ClC(Cl)Cl"

# Carbon tetrachloride; four Cl
# sml = "ClC(Cl)(Cl)Cl"

# sml = "c1ccccc1" # benzene; boring

# Calcium carbonate; boring
# sml = "[Ca+2].[O-]C([O-])=O"

# Potassium cyanide; three elements
# sml = "[K+].[C-]#N"

# Four atoms, each a different element
# sml = "C(=O)(F)N"

# FeCl3; good for distinguishing very similar masses
# sml = "Cl[Fe](Cl)Cl"

# Cl2; simple case
# sml = "ClCl"

# Homo triatomic
# sml = "[O-][O+]=O"

# Homo diatomic
# sml = "O=O"

# %%
mol = Chem.MolFromSmiles(sml)

# %%
mol

# %%
def composition(molecule):
    """Get the composition of an RDKit molecule:
    Atomic counts, including hydrogen atoms, and any charge.
    For example, fluoride ion (chemical formula F-, SMILES string [F-])
    returns {9: 1, 0: -1}.

    :param molecule: The molecule to analyze
    :type some_input: An RDKit molecule
    :rtype: A dictionary.
    """
    # Check that there is a valid molecule
    if molecule:

        # Add hydrogen atoms--RDKit excludes them by default
        molecule_with_Hs = Chem.AddHs(molecule)
        comp = defaultdict(lambda: 0)

        # Get atom counts
        for atom in molecule_with_Hs.GetAtoms():
            comp[atom.GetSymbol()] += 1

        return comp

# %%
def flatten_ndarray(arr):
    flat_list = []
    for element in arr:
        if isinstance(element, np.ndarray):
            flat_list.extend(flatten_ndarray(element))
        else:
            flat_list.append(element)
    return flat_list

# %%
def binomial_coefficient(n, k):
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# %%
def distribute_items(n, k):
    def backtrack(remaining_items, bins, current_bin):
        if current_bin == k:
            if remaining_items == 0:
                results.append(list(bins))
            return

        for items_in_current_bin in range(remaining_items + 1):
            bins[current_bin] = items_in_current_bin
            backtrack(remaining_items - items_in_current_bin, bins, current_bin + 1)

    results = []
    bins = [0] * k
    backtrack(n, bins, 0)
    return results

# %%
def atom_set(molecule):
    """Get the composition of an RDKit molecule:
    Atomic counts, including hydrogen atoms, and any charge.
    For example, fluoride ion (chemical formula F-, SMILES string [F-])
    returns {9: 1, 0: -1}.

    :param molecule: The molecule to analyze
    :type some_input: An RDKit molecule
    :rtype: A dictionary.
    """
    # Check that there is a valid molecule
    if molecule:
      # Add hydrogen atoms--RDKit excludes them by default
      molecule_with_Hs = Chem.AddHs(molecule)
      atom_symbols = set(atom.GetSymbol() for atom in molecule_with_Hs.GetAtoms())
      return atom_symbols

# %%
from mendeleev import element, isotope

set_of_atoms = atom_set(mol)
print(set_of_atoms)
isotopes = {}
for atom in set_of_atoms:
  # print(f"{atom=}")
  # print(f"{element(atom)=}")
  element_isotopes = element(atom).isotopes
  occuring_isotopes = [[isotope.mass_number, isotope.abundance/100] for isotope in element_isotopes if isotope.abundance != None]
  isotopes[atom] = occuring_isotopes
print(isotopes)

# %%
def print_ndarray_elements(arr, prefix=''):
    if isinstance(arr, Chem.Mol):
        print(prefix + Chem.MolToSmiles(arr))
    else:
        for i in range(len(arr)):
            new_prefix = prefix + "[" + str(i) + "]"
            print_ndarray_elements(arr[i], new_prefix)

# %%
def extract_properties(arr, prop_name):
    if isinstance(arr, np.ndarray):
        # If arr is a NumPy ndarray, iterate through its elements and apply the function recursively
        return np.array([extract_properties(item, prop_name) for item in arr])
    elif isinstance(arr, list):
        # If arr is a list, iterate through its elements and apply the function recursively
        return [extract_properties(item, prop_name) for item in arr]
    elif hasattr(arr, prop_name):
        # If arr has the specified property, extract its value
        return getattr(arr, prop_name)
    else:
        # If the property is not found, return None
        return None

# %%
def assign_isotopes(
      arr, 
      isotope_count_distribution, 
      prefix, 
      n_this_element,
      ):
    if isinstance(arr, np.ndarray):
        for i in range(len(arr)):
            new_prefix = prefix + [i]
            assign_isotopes(arr[i], isotope_count_distribution, new_prefix, n_this_element)
    else:
        # Clone the distribution, making it independent so the next loop won't be affected
        distribution = copy.deepcopy(isotope_count_distribution[prefix[0]])

        # Calculate abundances
        if arr.abundance is None:
           arr.abundance = 1
        a = 1
        b = binomial_coefficient(n_this_element, distribution[0])
        for isotope_index, isotope_count in enumerate(distribution):
          a_this_isotope = isotopes[this_element][isotope_index][1]
          a *= a_this_isotope**isotope_count
        arr.abundance *= a*b

        # print(f"{distribution=}")
        # Assign isotopes to atoms of this element type
        for atom_index, atom in enumerate(arr.mol.GetAtoms()):
          if atom.GetSymbol() == this_element:
            # Count down each isotope, going on to the next when zero
            for isotope_index, isotope in enumerate(distribution):
              # print(f"    {isotope_index=}")
              if distribution[isotope_index] > 0:
                  atom.SetIsotope(isotopes[this_element][isotope_index][0])
                  distribution[isotope_index] -= 1

        isotopes_list = []
        for atom_index, atom in enumerate(arr.mol.GetAtoms()):
          isotopes_list.append(atom.GetIsotope())
        isotopes_str = ','.join(map(str, isotopes_list))
        # print("".join([f"[{i}]" for i in prefix]) + Chem.MolToSmiles(arr.mol) + isotopes_str)
    return arr

# %%
def update_molecular_isotopes(
      arr, 
      ):
    if isinstance(arr, np.ndarray):
        for i in range(len(arr)):
            update_molecular_isotopes(arr[i])
    elif isinstance(arr, MolecularIsotope):
        arr.update()

# %%
# Initialize ndarray to hold MolecularIsotope objects
molecular_isotopes:np.ndarray = np.array(MolecularIsotope)

# Loop through the elements in this molecule
for this_element, n_this_element in composition(mol).items():
  n_isotopes_this_element = len(isotopes[this_element])

  # Calculate possible distributions of isotopes across atoms of this element
  isotope_count_distribution = distribute_items(n_this_element, n_isotopes_this_element)
  n_distributions = len(isotope_count_distribution)

  if molecular_isotopes.shape == ():
     # Start by creating a 1-dimensional array,
     # making sure to make each mol an independent object
     molecular_isotopes = np.array([MolecularIsotope(sml=sml, mol=Chem.Mol(mol)) for _ in range(n_distributions)], dtype=object)
  else:
    # Create a list of m copies of the current object, 
    # namely the n-1 dimensional ndarray representing elements prior to this element
    molecular_isotopes_list = [copy.deepcopy(molecular_isotopes) for _ in range(n_distributions)]
    
    # Convert the list of copies to a NumPy ndarray
    molecular_isotopes = np.array(molecular_isotopes_list, dtype=object)
    
  # Assign isotopes and abundances
  molecular_isotopes = assign_isotopes(
     arr=molecular_isotopes, 
     isotope_count_distribution=isotope_count_distribution, 
     prefix=[],
     n_this_element=n_this_element,
     )
  
# Update the properties of each MolecularIsotope to get exact masses
update_molecular_isotopes(molecular_isotopes)

# %%
print(f"{molecular_isotopes.shape=}")
mols_array = extract_properties(molecular_isotopes, "mol")
mols_flat = flatten_ndarray(mols_array)

abundances_array = extract_properties(molecular_isotopes, "abundance")
abundances_flat = flatten_ndarray(abundances_array)
masses_array = extract_properties(molecular_isotopes, "mass")
masses_flat = flatten_ndarray(masses_array)

abundances_flat_str = [f"{mass:.3f}: {abundance:.3f}" for mass, abundance in zip(masses_flat, abundances_flat)]

Chem.Draw.MolsToGridImage(mols_flat, legends=abundances_flat_str, subImgSize=(100, 100)) #, legends=[str(mass) for mass in masses])

# %% [markdown]
# Divider

# %%
# #Debugging only!
# # isotopes = {'Cl': [[35, 0.5], [37, 0.5]]}

# n_atoms = 2

# n_enumerated = 0
# # n_test_atoms = 4
# # for isotope_count in range(n_test_atoms + 1):
# for isotope_count in range(n_atoms + 1):
#     n_enumerated += binomial_coefficient(n_atoms, isotope_count)

# print(f"{n_enumerated=}")

# #Debugging only!
# # n_atoms = 3

# sum_ab = 0
# for isotope_count in range(n_atoms + 1):
#   # for isotope in isotopes["Cl"]:
#     n_this_isotope = isotope_count
#     mass_this_isotope = isotopes["Cl"][0][0]
#     n_other_isotope = n_atoms - isotope_count
#     mass_other_isotope = isotopes["Cl"][1][0]
#     a_this_isotope = isotopes["Cl"][0][1]
#     a_other_isotope = isotopes["Cl"][1][1]
#     b = binomial_coefficient(n_atoms, n_this_isotope)
#     a_from_this_isotope = (a_this_isotope**n_this_isotope)
#     a_from_other_isotope = (a_other_isotope**n_other_isotope)
#     # print(f"this isotope : {n_this_isotope} {mass_this_isotope} {a_this_isotope} {a_this_isotope**n_this_isotope} {binomial_coefficient(n_atoms, n_this_isotope)} {a_from_this_isotope}")
#     # print(f"other isotope: {n_other_isotope} {mass_other_isotope} {a_other_isotope} {a_other_isotope**n_other_isotope} {binomial_coefficient(n_atoms, n_other_isotope)} {a_from_other_isotope}")
#     headers = "n          Mass     A      Result       Frac A_from"

#     # Define the data for the first isotope
#     data1 = f"{n_this_isotope:<12} {mass_this_isotope:<6} {a_this_isotope:<6} {a_this_isotope**n_this_isotope:<12} {a_from_this_isotope}"

#     # Define the data for the other isotope
#     data2 = f"{n_other_isotope:<12} {mass_other_isotope:<6} {a_other_isotope:<6} {a_other_isotope**n_other_isotope:<12} {a_from_other_isotope}"

#     print(headers)
#     print(data1)
#     print(data2)

#     a = a_from_this_isotope * a_from_other_isotope
#     print(f"{a=} {b=} {a*b=}")
#     sum_ab += a*b
# print(f"{sum_ab=}")
#     # print(isotope)

# %%
sml_check = "C#C"
mol_check = Chem.MolFromSmiles(sml_check)
for atom_index, atom in enumerate(mol_check.GetAtoms()):
    print(f"{atom_index}, {atom.GetSymbol()}, {atom.GetIsotope()}")

# %%
for element, count in composition(mol).items():
    abundance = isotopic_abundances[element]
    contribution *= pow(abundance, count)
    contribution *= binomial_coefficient(molecular_formula[element_symbol], count)
return contribution



