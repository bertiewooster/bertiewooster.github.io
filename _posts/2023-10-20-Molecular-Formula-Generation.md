# Molecular Formula Generation

In cheminformatics, the typical way of representing a molecule is with a SMILES string such as `CCO` for ethanol. A SMILES string can be converted into a molecular graph, which can be used to determine molecular structure and related properties. However, there are still cases where the molecular formula such as C<sub>2</sub>H<sub>6</sub>O is useful. 

For example, the molecular formula is sufficient to determine the molecular mass, to calculate the predicted results from an elemental analysis, to get a sense for the elemental composition of a given molecule, or balancing chemical equations. One of the best-known [chemical equations describes photosynthesis](https://bio.libretexts.org/Bookshelves/Botany/Botany_(Ha_Morrow_and_Algiers)/04%3A_Plant_Physiology_and_Regulation/4.01%3A_Photosynthesis_and_Respiration/4.1.03%3A_Photosynthesis_Overview_and_Equation):

<img alt="Photosynthesis chemical equation: 6CO2 + 6H2O → C6H12O6 + 6O2" src="data:image/png;base64,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" />

I'm working on a blog post where I need to calculate the molecular formula including isotopes. I was unable to find a package that produced molecular formulas including isotopes from SMILES strings, so I wrote a function to do so. (I tried [chemparse](https://pypi.org/project/chemparse/) and did not succeed.)


```python
!pip install rdkit
```


```python
from collections import defaultdict
from IPython.display import display, Markdown

from rdkit import Chem
```

## Setting up the molecule

Because RDKit excludes hydrogen atoms by default, but we want to include hydrogens in our molecular formula, we tell RDKit to add hydrogens.


```python
# Ethanol with isotopes
sml = "[13CH3]C[18OH]"
mol = Chem.AddHs(Chem.MolFromSmiles(sml))
mol
```




    
![Ethanol molecule where the oxygen is the O-18 isotope and the carbon not bonded to the oxygen in is C-13 isotope](/images/2023-10-20-Molecular-Formula-Generation_files/2023-10-20-Molecular-Formula-Generation_8_0.png)
    



## Composition function

The RDKit does not include a composition function to give the number of each element (including isotopes if desired) in a molecule, so we make one here. It's based on an [algorithm suggested by @IchiruTake](https://github.com/rdkit/rdkit/discussions/5339).


```python
def composition(
        molecule: Chem.Mol,
        isotopes: bool = False,
        ) -> defaultdict:
    """Get the composition of an RDKit molecule:
    Atomic counts, including hydrogen atoms, and isotopes if requested.
    For example, ethanol (SMILES [13C](H)(H)(H)CO, formula C2H6O) returns:
      if isotopes = False (default): {'C': 2, 'O': 1, 'H': 6}.
      if isotopes = True: {'C': {13: 1, 0: 1}, 'O': {0: 1}, 'H': {0: 6}}.

    :param molecule: The RDKit molecule to analyze
    :param isotopes: Whether to include the isotope of each atom
    :returns: if isotopes = False (default): a dictionary of element:count entries;
              if isotopes = True: a nested dictionary of element:isotope:count entries.
    """
    # Check that there is a valid molecule
    if not molecule:
        return

    # Add hydrogen atoms--RDKit excludes them by default
    molecule = Chem.AddHs(molecule)
    comp = defaultdict(lambda: 0)

    # Get atom counts
    for atom in molecule.GetAtoms():
        element_symbol = atom.GetSymbol()
        # If isotopes not requested, simply count the number of atoms of each element
        if not isotopes:
            comp[element_symbol] += 1
        # If isotopes requested, count the number of each isotope of each element
        else:
            isotope = atom.GetIsotope()
            try:
                comp[element_symbol][isotope] += 1
            except:
                comp[element_symbol] = defaultdict(lambda: 0)
                comp[element_symbol][isotope] += 1
    return comp
```

With `isotopes=False`, we get the count of each element:


```python
composition(mol, isotopes=False)
```




    defaultdict(<function __main__.composition.<locals>.<lambda>()>,
                {'C': 2, 'O': 1, 'H': 6})



With `isotopes=True`, the first layer (dictionary) is elements, and the second layer is the isotope and the count:


```python
composition(mol, isotopes=True)
```




    defaultdict(<function __main__.composition.<locals>.<lambda>()>,
                {'C': defaultdict(<function __main__.composition.<locals>.<lambda>()>,
                             {13: 1, 0: 1}),
                 'O': defaultdict(<function __main__.composition.<locals>.<lambda>()>,
                             {18: 1}),
                 'H': defaultdict(<function __main__.composition.<locals>.<lambda>()>,
                             {0: 6})})



## Molecular formula generation

Now that we have the count of each element (with isotopes if desired) in the molecule, we can build the molecular formula. Because of the superscripts for isotopes, and subscripts for counts, we need a formatting language. Our `mol_to_formatted_formula()` provides a dictionary with two options:
- [Markdown](https://en.wikipedia.org/wiki/Markdown), which is commonly converted to HTML. Here the only formatting used are the HTML tags for superscripts and subscripts, so the outputs are valid HTML as well.
- [LaTeX](https://en.wikipedia.org/wiki/LaTeX), which is commonly used for technical papers. LaTeX would overlap a subscript and a superscript that immediately followed it, so we add a small amount of horizontal space using an empty group `{}`.


```python
def mol_to_formatted_formula(
        mol: Chem.Mol,
        isotopes: bool = False,
        ) -> dict[str, str]:
    if mol is not None:
        comp = composition(mol, isotopes)

        formula = {'markdown': "", 'latex': ""}

        if isotopes:
            isotopes_dict = defaultdict(lambda: defaultdict(str))
            subscripts = defaultdict(lambda: defaultdict(int))
            superscripts = defaultdict(list)

            for element, counts in comp.items():
                for isotope, count in counts.items():
                    if count > 1:
                        subscripts[element][isotope] = count
                    if isotope != 0:
                        superscripts[element].append(isotope)

                    isotopes_dict[element][isotope] = 1
                # Sort the element's isotopes from lowest to highest
                superscripts[element].sort()

            last_item_is_subscript = False
            sorted_element_keys = sorted(isotopes_dict.keys())
            for element in sorted_element_keys:
                isotope_count_pairs = isotopes_dict[element]
                # Sort the element's isotopes from lowest to highest
                sorted_isotope_keys = sorted(isotope_count_pairs.keys())
                for isotope in sorted_isotope_keys:
                    if element in superscripts:
                        if isotope in superscripts[element]:
                            formula["markdown"] += f"<sup>{isotope}</sup>"

                            # If superscript immediately follows subscript, 
                            # add a small amount of horizontal space using an empty group {}
                            # to prevent them from vertically overlapping
                            if last_item_is_subscript:
                                formula["latex"] += "{}"
                            formula["latex"] += "^{{{}}}".format(isotope)
                            last_item_is_subscript = False
                    formula["markdown"] += element
                    formula["latex"] += element
                    last_item_is_subscript = False
                    if element in subscripts:
                        if isotope in subscripts[element]:
                            isotope_count = subscripts[element][isotope]
                            formula["markdown"] += f"<sub>{isotope_count}</sub>"
                            formula["latex"] += "_{{{}}}".format(isotope_count)
                        last_item_is_subscript = True
            # Add beginning and ending dollar signs to LaTeX formula
            formula["latex"] = "$" + formula["latex"] + "$"
        else:
            # Handling the case when isotopes is False
            sorted_element_keys = sorted(comp.keys())

            for element in sorted_element_keys:
                count = comp[element]
                formula["markdown"] += element
                formula["latex"] += element
                if count > 1:
                    formula["markdown"] += f"<sub>{count}</sub>"
                    formula["latex"] += "_{{{}}}".format(count)
            formula["latex"] = "$" + formula["latex"] + "$"

        return formula
    else:
        return "Invalid molecule"
```

Here's the Markdown for our isotopic ethanol molecule:


```python
isotope_formula_markdown = mol_to_formatted_formula(mol, isotopes=True)["markdown"]
isotope_formula_markdown
```




    'C<sup>13</sup>CH<sub>6</sub><sup>18</sup>O'



and we display it as Markdown using `Markdown()`:


```python
Markdown(isotope_formula_markdown)
```




C<sup>13</sup>CH<sub>6</sub><sup>18</sup>O



Here's the LaTeX version:


```python
isotope_formula_latex = mol_to_formatted_formula(mol, isotopes=True)["latex"]
isotope_formula_latex
```




    '$C^{13}CH_{6}{}^{18}O$'



which can also be displayed using `Markdown()`:


```python
Markdown(isotope_formula_latex)
```




$C^{13}CH_{6}{}^{18}O$



As far as the order of elements in a chemical formula, `mol_to_formatted_formula()` simply alphabetizes them. So the elements some formulas may not appear in the typical order. Within an element, `mol_to_formatted_formula()` gives the isotopes in increasing order, with the unspecified isotope first.

To go directly from a SMILES string to a formula, we can use this utility function `smiles_to_formatted_formula()`:


```python
def smiles_to_formatted_formula(
        smiles:str,
        isotopes:bool=False):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol_to_formatted_formula(mol, isotopes=isotopes)
```


```python
isotope_formula_latex_from_smiles = smiles_to_formatted_formula(sml, isotopes=True)["latex"]
Markdown(isotope_formula_latex_from_smiles)
```




$C^{13}CH_{6}{}^{18}O$



### Improved formatting using LaTeX

LaTeX italicizes letters by default, so we can use the [LaTeX `\mathrm`](https://www.tutorialspoint.com/tex_commands/mathrm.htm) to make them non-italicized. (The C-style [string substitution using `%`](https://realpython.com/python-string-formatting/#1-old-style-string-formatting-operator) is necessary instead of f-string `{variable}` or `{}.format(variable)` because the braces in those last two interfere with the LaTeX formatting.)


```python
def markdown_formula(latex:str) -> str:
    latex_markdown = r"$\mathrm{ %s}$" % (latex.strip("$"))
    return latex_markdown
```

Here's the non-italicized result:


```python
Markdown(markdown_formula(isotope_formula_latex))
```




$\mathrm{ C^{13}CH_{6}{}^{18}O}$



As a further utility, we can immediately display the result as Markdown by incorporating that function:


```python
def display_markdown_formula(latex:str) -> str:
    latex_markdown = r"$\mathrm{ %s}$" % (latex.strip("$"))
    return Markdown(latex_markdown)
```


```python
display_markdown_formula(isotope_formula_latex)
```




$\mathrm{ C^{13}CH_{6}{}^{18}O}$



## Conclusion
Now that we have a way to calculate molecular formulas, and two formats to display them in, the next blog post will give an application of each format.

## Postscript
Here's how to use LaTeX to create the photosynthesis chemical equation shown at the top of the blog post.


```python
photosynthesis_smls = {
    "Water": "O",
    "Oxygen": "O=O",
    "Carbon Dioxide": "O=C=O",
    "Glucose": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
}
```


```python
mols = {name:Chem.MolFromSmiles(sml) for name, sml in photosynthesis_smls.items()}
```


```python
formulas = {name:markdown_formula(mol_to_formatted_formula(mol)["latex"]) for name, mol in mols.items()}
```


```python
photosynthesis = "$" + "6" + formulas["Carbon Dioxide"].strip("$") + "+ 6" + formulas["Water"].strip("$") + "→" + formulas["Glucose"].strip("$") + "+ 6" + formulas["Oxygen"].strip("$") + "$"
display(Markdown(photosynthesis))
```


$6\mathrm{CO_{2}}+ 6\mathrm{H_{2}O}→\mathrm{C_{6}H_{12}O_{6}}+ 6\mathrm{O_{2}}$

