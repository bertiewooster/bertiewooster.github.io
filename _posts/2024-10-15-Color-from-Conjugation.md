# Why some organic molecules have a color

## Correlating optical absorption wavelength with conjugated bond chain length

Molecules have a color if their electronic energy levels are close enough to absorb visible rather than ultraviolet light. For organic molecules, that's often because of an extensive chain of conjugated bonds. Can we use cheminformatics to find evidence that increasing conjugated bond chain length decreases absorption wavelength, which makes a molecule colored?

*[Open this notebook in Google Colab](https://colab.research.google.com/drive/19ZxvOVLYdCU5c47x8CZ1ns8MFyXTz9YV?usp=sharing) so you can run it without installing anything on your computer*

## Introduction

Colored molecules have applications in television screens, sensors, and to give color to fabrics, paints, foods, and more. I did my PhD in optical spectroscopy, so I was interested by the open-access database [Experimental database of optical properties of organic compounds](https://www.nature.com/articles/s41597-020-00634-8) with 20,236 data points. One of the first things that came to mind was that organic colored compounds often get their color from an extensive chain of conjugated bonds. Here are quotes from an [online textbook](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Spectroscopy/Electronic_Spectroscopy/Electronic_Spectroscopy_Basics/What_Causes_Molecules_to_Absorb_UV_and_Visible_Light):

> If you extend this to compounds with really massive delocalisation, the wavelength absorbed will eventually be high enough to be in the visible region of the spectrum, and the compound will then be seen as colored. A good example of this is the orange plant pigment, beta-carotene - present in carrots, for example.

> The more delocalization there is, the smaller the gap between the highest energy pi bonding orbital and the lowest energy pi anti-bonding orbital. To promote an electron therefore takes less energy in beta-carotene than in the cases we've looked at so far - because the gap between the levels is less.

Here are two examples of molecules, one cyclic and one acylic, with conjugated pi bonds. In each case, the p orbitals on adjacent atoms line up so that the electrons are delocalized over all the conjugated bonds in the molecule (which is all the bonds in these two molecules).

![Molecules and their conjugated pi bonds](/images/molecules_and_conjugated.png)

*Attribution: [Conjugated Pi Bond Systems from LibreTexts](https://chem.libretexts.org/Courses/Nassau_Community_College/Organic_Chemistry_I_and_II/02%3A_Structure_and_Properties_of_Organic_Molecules/2.04%3A_2.4_Conjugated_Pi_Bond_Systems), remixed by Jeremy Monat*

The more bonds that electrons can delocalize over, the more pi bonding and anti-bonding orbitals; and the highest occupied molecular orbital (HOMO, the top green line in each diagram below, the highest-energy pi bonding orbital) increases in energy, while the lowest unoccupied molecular orbital (LUMO, the lowest red line in each diagram below, the lowest-energy pi anti-bonding orbital) decreases in energy. The gap between the two becomes smaller, and if the conjugated chain is long enough, the HOMO-LUMO energy gap becomes small enough that it's in the visible spectrum rather than the ultraviolet. When a molecule absorbs visible light, we perceive that it has a color.

![Molecular orbital energy diagram for ethene, buta-1,3-diene, and hexa-1,3,5-triene](/images/conjugated_molecular_orbitals.gif)

*Attribution: [What Causes Molecules to Absorb UV and Visible Light from LibreTexts](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Spectroscopy/Electronic_Spectroscopy/Electronic_Spectroscopy_Basics/What_Causes_Molecules_to_Absorb_UV_and_Visible_Light), authored, remixed, and/or curated by [Jim Clark](http://www.chemguide.co.uk/).*

The visible spectrum starts at about 400 nm (violet) and goes to about 740 nm (red). So a molecule that absorbs light in that range will be perceived as colored.

![Visible spectrum](/images/Linear_visible_spectrum.svg)

*Attribution: [Gringer, Public domain, via Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Linear_visible_spectrum.svg).*

### Cheminformatics exploration

Can we use cheminformatics to find evidence that increasing conjugated bond chain length decreases absorption wavelength? To check, I used the open-access database [Experimental database of optical properties of organic compounds](https://www.nature.com/articles/s41597-020-00634-8) from 2020 with 20,236 data points. The [optical data can be downloaded as a CSV file](https://figshare.com/articles/dataset/DB_for_chromophore/12045567/2?file=23637518).

## Packages setup


```python
import math
from typing import Iterable, Dict, List, Tuple
from IPython.display import display, Math, SVG
from functools import cache
import warnings

from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Mol
from rdkit.Geometry import Point2D
import numpy as np
import polars as pl
import polars.selectors as cs
from polars.exceptions import ColumnNotFoundError
import altair as alt
from great_tables import GT, style, loc
import latexify
```


```python
# Suppress the specific RDKit IPythonConsole truncation warning for MolsToGridImage
warnings.filterwarnings(
    "ignore", message=r"Truncating the list of molecules to be displayed to \d+"
)
```

## Find longest conjugated bond chain

To find the longest conjugated bond chain in each molecule, we'll use a graph-traversal algorithm. We'll define two functions--one to find connected bonds and one to get the longest conjugated bond chain in a molecule--then describe how they work using an example molecule.


```python
def find_connected_bonds(
    graph: np.ndarray,
    start_bond: int,
    visited: Iterable,
) -> List:
    """
    Find connected bonds in an adjacency matrix. Terminology is specific to bonds in molecules, but algorithm should work for any adjacency matrix.

    :param graph: The bond adjacency matrix
    :param start_bond: The bond index to start the chain from
    :param visited: The bond indices already visited; preferably a set for efficiency
    :returns: The list of bonds connected to, and including, the start bond
    """
    # Initialize the stack with the start bond
    stack = [start_bond]
    connected_component = []

    while stack:
        # Pop off the stack the last bond added
        bond_index = stack.pop()
        # Only follow the chain from this bond if this bond hasn't already been visited
        if bond_index not in visited:
            # Add the bond index to the list of visited bonds
            visited.add(bond_index)
            # Note that the bond is connected to the start bond
            connected_component.append(bond_index)

            # Add all neighbors to the stack for traversal if they haven't already been visited
            for neighbor, is_connected in enumerate(graph[bond_index]):
                if is_connected and neighbor not in visited:
                    stack.append(neighbor)
    return connected_component
```


```python
def get_longest_conjugated_bond_chain(
    mol: Mol, verbose: bool = False
) -> tuple[List[int], np.ndarray]:
    """
    Get the longest conjugated bond chain in a molecule.

    :param mol: The RDKit molecule
    :param verbose: Whether to return the bond_matrix in addition to conjugated_bonds_out
    :returns: The list of bond indices of the longest conjugated bond chain in the molecule; if verbose is true, also the bond adjacency matrix
    """
    # Create a list to store conjugated bond indices
    conjugated_bonds = [
        bond.GetIdx() for bond in mol.GetBonds() if bond.GetIsConjugated()
    ]

    if not conjugated_bonds:
        # No conjugated bonds found, return empty list
        return []

    n_conjugated_bonds = len(conjugated_bonds)

    # Build a subgraph of the conjugated bonds only;
    #   initially populate it with zeroes, indicating bonds are not connected
    bond_matrix = np.zeros((n_conjugated_bonds, n_conjugated_bonds), dtype=int)

    # Populate the bond adjacency matrix for conjugated bonds
    for i, bond_i in enumerate(conjugated_bonds):
        bond_i_obj = mol.GetBondWithIdx(bond_i)
        for j, bond_j in enumerate(conjugated_bonds):
            if i != j:
                bond_j_obj = mol.GetBondWithIdx(bond_j)
                # Check if two conjugated bonds share an atom--
                #   do the set of {beginning atom, ending atom} overlap for the two bonds
                if (
                    len(
                        set([bond_i_obj.GetBeginAtomIdx(), bond_i_obj.GetEndAtomIdx()])
                        & set(
                            [bond_j_obj.GetBeginAtomIdx(), bond_j_obj.GetEndAtomIdx()]
                        )
                    )
                    > 0
                ):
                    # Change the bond matrix value to 1, indicating the two bonds are connected
                    bond_matrix[i, j] = 1
                    bond_matrix[j, i] = 1

    # Initialize variables to store the longest conjugated bond chain
    visited = set()
    longest_bond_chain = []

    # Starting from each bond, traverse the graph and find the largest connected component
    for start_bond in range(n_conjugated_bonds):
        if start_bond not in visited:
            bond_chain = find_connected_bonds(bond_matrix, start_bond, visited)
            # Note that bonds are added to `visited` in find_connected_bonds(),
            #   so any bonds already visited from a previous starting bond
            #   won't have find_connected_bonds run on it
            # If this chain is longer than the longest one found so far, mark this chain as the longest
            if len(bond_chain) > len(longest_bond_chain):
                longest_bond_chain = bond_chain

    # Convert subgraph bond indices back to the original bond indices
    conjugated_bonds_out = [conjugated_bonds[i] for i in longest_bond_chain]
    conjugated_bonds_out.sort()

    if not verbose:
        return conjugated_bonds_out
    else:
        return conjugated_bonds_out, bond_matrix
```

Let's use an example branched molecule to explain how these algorithms work.


```python
C6H8 = Chem.MolFromSmiles("C=CC(=C)C=C")
longest_conjugated_bond_chain, bond_matrix = get_longest_conjugated_bond_chain(
    mol=C6H8, verbose=True
)
Draw.MolToImage(
    C6H8,
    highlightBonds=longest_conjugated_bond_chain,
)
```




    
![Branched compound 3-methylidenepenta-1,4-diene with all bonds highlighted because they are all conjugated](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_25_0.png)
    



Let's label the bond indices using the following function.


```python
def label_bonds(mol: Mol, offset_y: float = 0) -> str:
    """
    Label the bonds with their indices in a molecule.

    :param mol: The RDKit molecule
    :param offset_y: The vertical offset for the bond labels; helpful for large molecules where the bond index numbers can be difficult to read
    :returns: SVG of the labeled molecule; display in a Jupyter Notebook using display(SVG()), e.g. display(SVG(label_bonds(my_mol)))
    """
    # Generate 2D coordinates for visualization
    Chem.rdDepictor.Compute2DCoords(mol)

    # Define image size and initialize an SVG drawer
    drawer = Draw.MolDraw2DSVG(300, 150)

    # Draw the molecule first
    drawer.DrawMolecule(mol)

    # Add bond numbers
    for bond in mol.GetBonds():
        idx = bond.GetIdx()

        # Get the positions of the atoms of the bond
        begin_atom_pos = mol.GetConformer().GetAtomPosition(bond.GetBeginAtomIdx())
        end_atom_pos = mol.GetConformer().GetAtomPosition(bond.GetEndAtomIdx())

        # Calculate midpoint of bond as midpoint between the two atoms
        mid_x = (begin_atom_pos.x + end_atom_pos.x) / 2
        mid_y = (begin_atom_pos.y + end_atom_pos.y) / 2

        # Optionally, offset the y coordinate to move the label just above (positive offset) or below (negative offset) the bond
        mid_point = Point2D(mid_x, mid_y + offset_y)

        # Add bond index at the offset midpoint
        drawer.DrawString(str(idx), mid_point)

    drawer.FinishDrawing()

    # Get the SVG text
    svg = drawer.GetDrawingText()

    return svg
```


```python
C6H8_img = label_bonds(C6H8)
display(SVG(C6H8_img))
```


    
![Branched compound 3-methylidenepenta-1,4-diene with the bonds numbered: 0, 1, 3, 4 along the molecule's spine, and 2 for the branch between bonds 1 and 3](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_28_0.svg)
    


To understand the bond adjacency matrix, let's use [Polars' Great Tables integration](https://posit-dev.github.io/great-tables/blog/polars-styling/) to make a nicely-formatted table.


```python
# Convert the bond matrix to a Polars DataFrame
df_bond_matrix = pl.DataFrame(bond_matrix)

# Rename columns to add custom labels
df_bond_matrix = df_bond_matrix.rename(
    {f"column_{str(i)}": f"Bond {i}" for i in range(bond_matrix.shape[1])}
)

# Add row index labels (Polars doesn't have an index column, so we'll add it as a new column)
row_indices = [f"Bond {i}" for i in range(bond_matrix.shape[0])]
index_col = "Adjacent?"
df_bond_matrix = df_bond_matrix.insert_column(0, pl.Series(index_col, row_indices))

# Add a Total column to sum up how many bonds a given bond is connected to
total_col = "Total"
df_bond_matrix = df_bond_matrix.with_columns(
    pl.sum_horizontal(cs.starts_with("Bond")).alias(total_col)
)

# Use GreatTables to format the table
GT(df_bond_matrix).tab_options(
    # Bold the column headings
    column_labels_font_weight="bold",
).tab_style(
    # Bold the index and total columns
    style=style.text(weight="bold"),
    locations=loc.body(columns=[index_col, total_col]),
).tab_header(
    # Add a title to the table
    title="Bond adjacency matrix"
)
```




<div id="dminubznkh" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>
#dminubznkh table {
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
          -webkit-font-smoothing: antialiased;
          -moz-osx-font-smoothing: grayscale;
        }

#dminubznkh thead, tbody, tfoot, tr, td, th { border-style: none !important; }
 tr { background-color: transparent !important; }
#dminubznkh p { margin: 0 !important; padding: 0 !important; }
 #dminubznkh .gt_table { display: table !important; border-collapse: collapse !important; line-height: normal !important; margin-left: auto !important; margin-right: auto !important; color: #333333 !important; font-size: 16px !important; font-weight: normal !important; font-style: normal !important; background-color: #FFFFFF !important; width: auto !important; border-top-style: solid !important; border-top-width: 2px !important; border-top-color: #A8A8A8 !important; border-right-style: none !important; border-right-width: 2px !important; border-right-color: #D3D3D3 !important; border-bottom-style: solid !important; border-bottom-width: 2px !important; border-bottom-color: #A8A8A8 !important; border-left-style: none !important; border-left-width: 2px !important; border-left-color: #D3D3D3 !important; }
 #dminubznkh .gt_caption { padding-top: 4px !important; padding-bottom: 4px !important; }
 #dminubznkh .gt_title { color: #333333 !important; font-size: 125% !important; font-weight: initial !important; padding-top: 4px !important; padding-bottom: 4px !important; padding-left: 5px !important; padding-right: 5px !important; border-bottom-color: #FFFFFF !important; border-bottom-width: 0 !important; }
 #dminubznkh .gt_subtitle { color: #333333 !important; font-size: 85% !important; font-weight: initial !important; padding-top: 3px !important; padding-bottom: 5px !important; padding-left: 5px !important; padding-right: 5px !important; border-top-color: #FFFFFF !important; border-top-width: 0 !important; }
 #dminubznkh .gt_heading { background-color: #FFFFFF !important; text-align: center !important; border-bottom-color: #FFFFFF !important; border-left-style: none !important; border-left-width: 1px !important; border-left-color: #D3D3D3 !important; border-right-style: none !important; border-right-width: 1px !important; border-right-color: #D3D3D3 !important; }
 #dminubznkh .gt_bottom_border { border-bottom-style: solid !important; border-bottom-width: 2px !important; border-bottom-color: #D3D3D3 !important; }
 #dminubznkh .gt_col_headings { border-top-style: solid !important; border-top-width: 2px !important; border-top-color: #D3D3D3 !important; border-bottom-style: solid !important; border-bottom-width: 2px !important; border-bottom-color: #D3D3D3 !important; border-left-style: none !important; border-left-width: 1px !important; border-left-color: #D3D3D3 !important; border-right-style: none !important; border-right-width: 1px !important; border-right-color: #D3D3D3 !important; }
 #dminubznkh .gt_col_heading { color: #333333 !important; background-color: #FFFFFF !important; font-size: 100% !important; font-weight: bold !important; text-transform: inherit !important; border-left-style: none !important; border-left-width: 1px !important; border-left-color: #D3D3D3 !important; border-right-style: none !important; border-right-width: 1px !important; border-right-color: #D3D3D3 !important; vertical-align: bottom !important; padding-top: 5px !important; padding-bottom: 5px !important; padding-left: 5px !important; padding-right: 5px !important; overflow-x: hidden !important; }
 #dminubznkh .gt_column_spanner_outer { color: #333333 !important; background-color: #FFFFFF !important; font-size: 100% !important; font-weight: bold !important; text-transform: inherit !important; padding-top: 0 !important; padding-bottom: 0 !important; padding-left: 4px !important; padding-right: 4px !important; }
 #dminubznkh .gt_column_spanner_outer:first-child { padding-left: 0 !important; }
 #dminubznkh .gt_column_spanner_outer:last-child { padding-right: 0 !important; }
 #dminubznkh .gt_column_spanner { border-bottom-style: solid !important; border-bottom-width: 2px !important; border-bottom-color: #D3D3D3 !important; vertical-align: bottom !important; padding-top: 5px !important; padding-bottom: 5px !important; overflow-x: hidden !important; display: inline-block !important; width: 100% !important; }
 #dminubznkh .gt_spanner_row { border-bottom-style: hidden !important; }
 #dminubznkh .gt_group_heading { padding-top: 8px !important; padding-bottom: 8px !important; padding-left: 5px !important; padding-right: 5px !important; color: #333333 !important; background-color: #FFFFFF !important; font-size: 100% !important; font-weight: initial !important; text-transform: inherit !important; border-top-style: solid !important; border-top-width: 2px !important; border-top-color: #D3D3D3 !important; border-bottom-style: solid !important; border-bottom-width: 2px !important; border-bottom-color: #D3D3D3 !important; border-left-style: none !important; border-left-width: 1px !important; border-left-color: #D3D3D3 !important; border-right-style: none !important; border-right-width: 1px !important; border-right-color: #D3D3D3 !important; vertical-align: middle !important; text-align: left !important; }
 #dminubznkh .gt_empty_group_heading { padding: 0.5px !important; color: #333333 !important; background-color: #FFFFFF !important; font-size: 100% !important; font-weight: initial !important; border-top-style: solid !important; border-top-width: 2px !important; border-top-color: #D3D3D3 !important; border-bottom-style: solid !important; border-bottom-width: 2px !important; border-bottom-color: #D3D3D3 !important; vertical-align: middle !important; }
 #dminubznkh .gt_from_md> :first-child { margin-top: 0 !important; }
 #dminubznkh .gt_from_md> :last-child { margin-bottom: 0 !important; }
 #dminubznkh .gt_row { padding-top: 8px !important; padding-bottom: 8px !important; padding-left: 5px !important; padding-right: 5px !important; margin: 10px !important; border-top-style: solid !important; border-top-width: 1px !important; border-top-color: #D3D3D3 !important; border-left-style: none !important; border-left-width: 1px !important; border-left-color: #D3D3D3 !important; border-right-style: none !important; border-right-width: 1px !important; border-right-color: #D3D3D3 !important; vertical-align: middle !important; overflow-x: hidden !important; }
 #dminubznkh .gt_stub { color: #333333 !important; background-color: #FFFFFF !important; font-size: 100% !important; font-weight: initial !important; text-transform: inherit !important; border-right-style: solid !important; border-right-width: 2px !important; border-right-color: #D3D3D3 !important; padding-left: 5px !important; padding-right: 5px !important; }
 #dminubznkh .gt_stub_row_group { color: #333333 !important; background-color: #FFFFFF !important; font-size: 100% !important; font-weight: initial !important; text-transform: inherit !important; border-right-style: solid !important; border-right-width: 2px !important; border-right-color: #D3D3D3 !important; padding-left: 5px !important; padding-right: 5px !important; vertical-align: top !important; }
 #dminubznkh .gt_row_group_first td { border-top-width: 2px !important; }
 #dminubznkh .gt_row_group_first th { border-top-width: 2px !important; }
 #dminubznkh .gt_striped { background-color: rgba(128,128,128,0.05) !important; }
 #dminubznkh .gt_table_body { border-top-style: solid !important; border-top-width: 2px !important; border-top-color: #D3D3D3 !important; border-bottom-style: solid !important; border-bottom-width: 2px !important; border-bottom-color: #D3D3D3 !important; }
 #dminubznkh .gt_sourcenotes { color: #333333 !important; background-color: #FFFFFF !important; border-bottom-style: none !important; border-bottom-width: 2px !important; border-bottom-color: #D3D3D3 !important; border-left-style: none !important; border-left-width: 2px !important; border-left-color: #D3D3D3 !important; border-right-style: none !important; border-right-width: 2px !important; border-right-color: #D3D3D3 !important; }
 #dminubznkh .gt_sourcenote { font-size: 90% !important; padding-top: 4px !important; padding-bottom: 4px !important; padding-left: 5px !important; padding-right: 5px !important; text-align: left !important; }
 #dminubznkh .gt_left { text-align: left !important; }
 #dminubznkh .gt_center { text-align: center !important; }
 #dminubznkh .gt_right { text-align: right !important; font-variant-numeric: tabular-nums !important; }
 #dminubznkh .gt_font_normal { font-weight: normal !important; }
 #dminubznkh .gt_font_bold { font-weight: bold !important; }
 #dminubznkh .gt_font_italic { font-style: italic !important; }
 #dminubznkh .gt_super { font-size: 65% !important; }
 #dminubznkh .gt_footnote_marks { font-size: 75% !important; vertical-align: 0.4em !important; position: initial !important; }
 #dminubznkh .gt_asterisk { font-size: 100% !important; vertical-align: 0 !important; }

</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
<thead>

  <tr class="gt_heading">
    <td colspan="7" class="gt_heading gt_title gt_font_normal">Bond adjacency matrix</td>
  </tr>
<tr class="gt_col_headings">
  <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Adjacent?">Adjacent?</th>
  <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Bond 0">Bond 0</th>
  <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Bond 1">Bond 1</th>
  <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Bond 2">Bond 2</th>
  <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Bond 3">Bond 3</th>
  <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Bond 4">Bond 4</th>
  <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Total">Total</th>
</tr>
</thead>
<tbody class="gt_table_body">
  <tr>
    <td style="font-weight: bold;" class="gt_row gt_left">Bond 0</td>
    <td class="gt_row gt_right">0</td>
    <td class="gt_row gt_right">1</td>
    <td class="gt_row gt_right">0</td>
    <td class="gt_row gt_right">0</td>
    <td class="gt_row gt_right">0</td>
    <td style="font-weight: bold;" class="gt_row gt_right">1</td>
  </tr>
  <tr>
    <td style="font-weight: bold;" class="gt_row gt_left">Bond 1</td>
    <td class="gt_row gt_right">1</td>
    <td class="gt_row gt_right">0</td>
    <td class="gt_row gt_right">1</td>
    <td class="gt_row gt_right">1</td>
    <td class="gt_row gt_right">0</td>
    <td style="font-weight: bold;" class="gt_row gt_right">3</td>
  </tr>
  <tr>
    <td style="font-weight: bold;" class="gt_row gt_left">Bond 2</td>
    <td class="gt_row gt_right">0</td>
    <td class="gt_row gt_right">1</td>
    <td class="gt_row gt_right">0</td>
    <td class="gt_row gt_right">1</td>
    <td class="gt_row gt_right">0</td>
    <td style="font-weight: bold;" class="gt_row gt_right">2</td>
  </tr>
  <tr>
    <td style="font-weight: bold;" class="gt_row gt_left">Bond 3</td>
    <td class="gt_row gt_right">0</td>
    <td class="gt_row gt_right">1</td>
    <td class="gt_row gt_right">1</td>
    <td class="gt_row gt_right">0</td>
    <td class="gt_row gt_right">1</td>
    <td style="font-weight: bold;" class="gt_row gt_right">3</td>
  </tr>
  <tr>
    <td style="font-weight: bold;" class="gt_row gt_left">Bond 4</td>
    <td class="gt_row gt_right">0</td>
    <td class="gt_row gt_right">0</td>
    <td class="gt_row gt_right">0</td>
    <td class="gt_row gt_right">1</td>
    <td class="gt_row gt_right">0</td>
    <td style="font-weight: bold;" class="gt_row gt_right">1</td>
  </tr>
</tbody>


</table>

</div>




The table demonstrates which bonds are adjacent, in the sense that the two bonds share an atom. For example, bond 1 is adjacent to bonds 0, 2, and 3. That makes sense based on the molecular diagram.

The bond chain starts with a start bond, in this case 0, and follows all its adjacent bonds to make a chain. Here, the algorithm went to bond 1 (the only bond connected to bond 0), then at the branch chose to go off the molecule's spine (longest atom chain) to go to bond 2, then followed the other branch to complete the molecule's spine (bonds 3 and 4):


```python
longest_conjugated_bond_chain
```




    [0, 1, 2, 3, 4]



In this case, a single bond chain found all the conjugated bonds in the molecule. The algorithm loops over all conjugated bonds to make sure it finds the longest chain. But if a given start bond (e.g., 1) was already visited because it was in the bond chain of a previous starting bond (e.g., 0), the algorithm doesn't re-trace the same bond chain. This is a big computational savings because it avoids unnecessary graph traversals, which are expensive. This is facilitated by the variable `visited` being passed from `get_longest_conjugated_bond_chain()` to `find_connected_bonds()`, where it is modified by adding the nodes visited.

Additional computational savings comes from excluding the non-conjugated bonds from the adjacency matrix. While not noticeable for the example molecule because all its bonds are conjugated, this can greatly reduce the adjacency matrix, and thus the graph traversal, for molecules where some bonds are not conjugated.

### Check that longest conjugated bond chain gives expected results

Let's make sure that our algorithm gives the expected results for a variety of cyclic and acyclic molecules.


```python
examples = {
    "benzene": "c1ccccc1",
    "naphthalene": "c1c2ccccc2ccc1",
    "anthracene": "c1ccc2cc3ccccc3cc2c1",
    "toluene": "c1ccccc1C",
    "benzene + 5 linear conjugated bonds": "c1ccccc1CC=CC=CC=C",
    "benzene + 7 linear conjugated bonds": "c1ccccc1CC=CC=CC=CC=C",
    "1,3-butadiene": "C=CC=C",
    "branched": "C=C\C=C/C(/C=C)=C/C=C",
    "branched + distant": "C=C\C=C\C\C=C\C=C/C(/C=C)=C/C=C",
}
mols = [Chem.MolFromSmiles(sml) for sml in examples.values()]
conjugated_bonds = [get_longest_conjugated_bond_chain(mol) for mol in mols]
Draw.MolsToGridImage(
    mols=mols,
    legends=examples.keys(),
    highlightBondLists=conjugated_bonds,
    subImgSize=(300, 200),
)
```




    
![Grid of molecules with each's longest conjugated chain highlighted](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_38_0.png)
    



Those highlighted conjugated bond chains are as expected, for example
- all the bonds are conjugated in benzene, as well as fused polycyclic molecules naphthalene and anthracene
- the bond off the ring in toluene is not conjugated
- if we put a side-chain on benzene with fewer than six conjugated bonds, the benzyl moiety remains the longest conjugated chain; but if we have seven conjugated bonds on the side chain, it becomes the longest conjugated chain
- in the "branched + distant" molecule, if we break the conjugation chain by having two C-C single bonds in a row, the chain does not include the distant, disconnected conjugated bonds

Of course we should check [beta-carotene](https://en.wikipedia.org/wiki/%CE%92-Carotene), whose color comes from its extended conjugated bond chain and is "[responsible for the orange color of carrots](https://chem.libretexts.org/Courses/Nassau_Community_College/Organic_Chemistry_I_and_II/02%3A_Structure_and_Properties_of_Organic_Molecules/2.04%3A_2.4_Conjugated_Pi_Bond_Systems)". beta-carotene "[absorbs most strongly between 400-500 nm. This is the green/blue part of the spectrum.](https://www.chm.bris.ac.uk/motm/carotene/beta-carotene_colourings.html)" Because it absorbs those wavelengths, when we look at a sample of beta-carotene we see the reflected light of the complimentary colors, so we perceive an orange color. Note that beta-carotene's absorption isn't that much lower energy (higher wavelength) than ultraviolet (which goes up to 400 nm), so molecules with longer conjugated chains might absorb at higher wavelengths such as 565-590 nm (yellow), giving them perceived colors towards blue, magenta, and purple.


```python
beta_carotene = Chem.MolFromSmiles(
    "CC2(C)CCCC(\C)=C2\C=C\C(\C)=C\C=C\C(\C)=C\C=C\C=C(/C)\C=C\C=C(/C)\C=C\C1=C(/C)CCCC1(C)C"
)
longest_conjugated_bond_chain = get_longest_conjugated_bond_chain(beta_carotene)
print(
    f"beta-carotene's longest conjugated bond chain is {len(longest_conjugated_bond_chain)} bonds long."
)
Draw.MolToImage(
    beta_carotene,
    highlightBonds=longest_conjugated_bond_chain,
    size=(800, 200),
)
```

    beta-carotene's longest conjugated bond chain is 21 bonds long.





    
![beta-carotene molecular structure with its conjugated bond chain highlighted](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_41_1.png)
    



And beta-carotene's molecular structure is also beautifully symmetric.

Now that we have an algorithm to find the longest conjugated bond chain in a molecule, let's apply it to the optical dataset.

## Prepare data

Let's read into a Polars dataframe the data from the [optical dataset CSV file](<https://figshare.com/articles/dataset/DB_for_chromophore/12045567/2?file=23637518).>).


```python
df = pl.read_csv("../data/DB for chromophore_Sci_Data_rev02.csv")
df
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (20_236, 14)</small><table border="1" class="dataframe"><thead><tr><th>Tag</th><th>Chromophore</th><th>Solvent</th><th>Absorption max (nm)</th><th>Emission max (nm)</th><th>Lifetime (ns)</th><th>Quantum yield</th><th>log(e/mol-1 dm3 cm-1)</th><th>abs FWHM (cm-1)</th><th>emi FWHM (cm-1)</th><th>abs FWHM (nm)</th><th>emi FWHM (nm)</th><th>Molecular weight (g mol-1)</th><th>Reference</th></tr><tr><td>i64</td><td>str</td><td>str</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>str</td></tr></thead><tbody><tr><td>1</td><td>&quot;N#Cc1cc2ccc(O)cc2oc1=O&quot;</td><td>&quot;O&quot;</td><td>355.0</td><td>410.0</td><td>2.804262</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>187.1537</td><td>&quot;DOI: 10.1021/acs.jpcb.5b09905&quot;</td></tr><tr><td>2</td><td>&quot;N#Cc1cc2ccc([O-])cc2oc1=O&quot;</td><td>&quot;O&quot;</td><td>408.0</td><td>450.0</td><td>3.961965</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>43.0</td><td>186.14576</td><td>&quot;DOI: 10.1021/acs.jpcb.5b09905&quot;</td></tr><tr><td>3</td><td>&quot;CCCCCCCCCCCC#CC#CCCCCCCCCCN1C(…</td><td>&quot;ClC(Cl)Cl&quot;</td><td>526.0</td><td>535.0</td><td>3.602954</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1061.54348</td><td>&quot;https://doi.org/10.1002/smll.2…</td></tr><tr><td>4</td><td>&quot;[O-]c1c(-c2nc3ccccc3s2)cc2ccc3…</td><td>&quot;CC#N&quot;</td><td>514.0</td><td>553.72</td><td>3.81</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>67.4</td><td>350.42028</td><td>&quot;https://doi.org/10.1016/j.snb.…</td></tr><tr><td>5</td><td>&quot;[O-]c1c(-c2nc3ccccc3s2)cc2ccc3…</td><td>&quot;CS(C)=O&quot;</td><td>524.0</td><td>555.0</td><td>4.7</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>58.0</td><td>50.0</td><td>350.42028</td><td>&quot;https://doi.org/10.1016/j.snb.…</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>20232</td><td>&quot;N#Cc1c(N2CCCC2)cc(-c2ccccc2)c2…</td><td>&quot;C1CCOC1&quot;</td><td>358.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>350.41992</td><td>&quot;DOI: 10.1021/ol9000679&quot;</td></tr><tr><td>20233</td><td>&quot;N#Cc1c(N2CCCCC2)cc(-c2ccccc2)c…</td><td>&quot;C1CCOC1&quot;</td><td>348.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>364.4468</td><td>&quot;DOI: 10.1021/ol9000679&quot;</td></tr><tr><td>20234</td><td>&quot;N#Cc1c(N2CCCCC2)cc(-c2cccc3ccc…</td><td>&quot;C1CCOC1&quot;</td><td>344.0</td><td>460.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>414.50668</td><td>&quot;DOI: 10.1021/ol9000679&quot;</td></tr><tr><td>20235</td><td>&quot;N#Cc1c(N2CCCCC2)cc(-c2ccc3cccc…</td><td>&quot;C1CCOC1&quot;</td><td>346.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>414.50668</td><td>&quot;DOI: 10.1021/ol9000679&quot;</td></tr><tr><td>20236</td><td>&quot;N#Cc1c(N2CCCCC2)cc(-c2ccc3ccc4…</td><td>&quot;C1CCOC1&quot;</td><td>344.0</td><td>473.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>488.58856</td><td>&quot;DOI: 10.1021/ol9000679&quot;</td></tr></tbody></table></div>



The dataset gives absorption and emission maxima as wavelengths. That makes sense because spectroscopists wavelength describes the color of the light used in the laboratory. But to compare, for example, the absorption and emission maximum for a molecule, it's better to use energy units to express the difference between different molecular energy levels. So let's convert wavelengths to energies in [electron volts](https://en.wikipedia.org/wiki/Electronvolt), eV.

To determine the conversion factor, let's use the equation for energy `E` as a function of wavelength `λ`.


```python
@latexify.function
def E(
    h: float,
    c: float,
    λ: float,
) -> float:
    """Calculate the binomial coefficient: how many ways there are to choose k items from n items.

    :param h: Planck's constant
    :param c: speed of light
    :param λ: wavelength
    :returns: energy
    """
    return h * c / λ


E
```




$$ \displaystyle E(h, c, λ) = \frac{h c}{λ} $$



Let's plug in the values of the [physical constants](https://en.wikipedia.org/wiki/Physical_constant#Table_of_physical_constants) to four decimal places and use factor-label dimensional analysis:


```python
eqn = r"$E = \frac{6.6261 \times 10^{-34} \, \text{J} \cdot \text{s} \cdot 2.9979 \times 10^8 \, \text{m/s}}{\lambda \times 10^{-9} \, \text{m}} \cdot \frac{1 \, \text{eV}}{1.6022 \times 10^{-19} \, \text{J}} = \frac{1239.8 \, \text{eV}}{\lambda}$"
display(Math(eqn))
```


$\displaystyle E = \frac{6.6261 \times 10^{-34} \, \text{J} \cdot \text{s} \cdot 2.9979 \times 10^8 \, \text{m/s}}{\lambda \times 10^{-9} \, \text{m}} \cdot \frac{1 \, \text{eV}}{1.6022 \times 10^{-19} \, \text{J}} = \frac{1239.8 \, \text{eV}}{\lambda}$


Doing that in Python to store the value in a variable:


```python
h = 6.6261e-34  # J*s
c = 2.9979e8  # m/s
nm = 1e-9  # m
eV = 1.6022e-19  # J
eV_nm = h * c / (nm * eV)
eV_nm
```




    1239.8193228061414



Now we can convert absorption and emission maxima to energy units of eV, then calculate their difference as the [Stokes shift](https://en.wikipedia.org/wiki/Stokes_shift). The Stokes shift reflects how much the molecule relaxes from its initial excited state (the [Franck–Condon state](https://en.wikipedia.org/wiki/Franck%E2%80%93Condon_principle)) to the lowest vibrational level in the excited state that it typically emits light from. In the following diagram, the blue arrow represents absorption from the ground state to the Franck-Condon state, and the green arrow represents emission from the relaxed excited state back to the ground state. The blue arrow is longer, representing the greater energy of absorption. The difference between the vertical length of the blue and green arrows is the Stokes shift.

![Franck-Condon state and relaxation](/images/Franck_Condon_Diagram.svg)

*Attribution: [Franck Condon Diagram on Wikipedia](https://en.wikipedia.org/wiki/Franck%E2%80%93Condon_principle#/media/File:Franck_Condon_Diagram.svg) by Samoza, licensed under the Creative Commons Attribution-Share Alike 3.0 Unported license*


```python
# To prevent duplicate-column errors when re-running this code, drop the columns we're about to add
for column in [
    "longest_bond_indices",
    "Longest conjugated bond length",
    "Absorption max (eV)",
    "Emission max (eV)",
    "Stokes shift (eV)",
]:
    try:
        df.drop(column)
    except ColumnNotFoundError:
        pass
```

Now we can calculate the energies in eV from the wavelengths in nm.


```python
df = df.with_columns(
    [
        (eV_nm / pl.col("Absorption max (nm)")).alias("Absorption max (eV)"),
        (eV_nm / pl.col("Emission max (nm)")).alias("Emission max (eV)"),
    ]
).with_columns(
    (pl.col("Absorption max (eV)") - pl.col("Emission max (eV)")).alias(
        "Stokes shift (eV)"
    ),
)
```

## Finding the longest conjugated chain for each molecule

Now we come to the computationally-intensive operation: Finding the the longest conjugated chain for each molecule using `conjugated_chain()` which calls `get_longest_conjugated_bond_chain()`. We define a function that will return the indices and length of the longest conjugated bond chain. The data set repeats some chromophores: there are a little less than three rows per chromophore. So we'll cache the results for each chromophore to avoid recalculating for each of its rows. Python's built-in module `functools` include a [`cache`](https://docs.python.org/3/library/functools.html#functools.cache) decorator that makes this simple.


```python
@cache
def conjugated_chain(sml) -> Dict[str : List[int], str:int]:
    """
    Find the indices and length for the longest bond chain in a SMILES.

    :param sml: SMILES to be made into a molecule
    :returns: A dictionary of longest bond chain indices and longest bond chain length
    """
    return_dict = dict()
    mol = Chem.MolFromSmiles(sml)
    longest_bond_indices = get_longest_conjugated_bond_chain(mol)
    return_dict["longest_bond_indices"] = longest_bond_indices
    return_dict["Longest conjugated bond length"] = len(longest_bond_indices)
    return return_dict
```

Now we use [Polars' `map_elements`](https://docs.pola.rs/api/python/stable/reference/expressions/api/polars.Expr.map_elements.html) to calculate the longest bond chain for each molecule. Because `conjugated_chain()` returns a dictionary, Polars treats it as a [`struct`](https://docs.pola.rs/api/python/stable/reference/series/struct.html), which we can then [`unnest`](https://docs.pola.rs/api/python/stable/reference/series/api/polars.Series.struct.unnest.html) to create a column for each dictionary key-value pair. We then sort the dataframe to put the longest bond chain lengths first so we can examine those molecules. Finally, we use Polars' [`shrink_to_fit()`](https://docs.pola.rs/api/python/stable/reference/dataframe/api/polars.DataFrame.shrink_to_fit.html) to decrease the dataframe memory usage and prevent problems with plotting.

This entire operation takes about 12 seconds on my laptop with caching, and about 31 seconds without, which roughly reflects the ratio of unique chromophores to data points. Caching is thus demonstrated to be effective here.


```python
# This may take 12-60 seconds: Finding the longest conjugated chain in each molecule in the dataframe
df = (
    df.with_columns(
        conjugated=pl.col("Chromophore").map_elements(
            lambda sml: conjugated_chain(sml), return_dtype=pl.Struct
        )
    )
    .unnest("conjugated")
    .sort(pl.col("Longest conjugated bond length"), descending=True)
    .shrink_to_fit()
)
```


```python
df.head()
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (5, 19)</small><table border="1" class="dataframe"><thead><tr><th>Tag</th><th>Chromophore</th><th>Solvent</th><th>Absorption max (nm)</th><th>Emission max (nm)</th><th>Lifetime (ns)</th><th>Quantum yield</th><th>log(e/mol-1 dm3 cm-1)</th><th>abs FWHM (cm-1)</th><th>emi FWHM (cm-1)</th><th>abs FWHM (nm)</th><th>emi FWHM (nm)</th><th>Molecular weight (g mol-1)</th><th>Reference</th><th>Absorption max (eV)</th><th>Emission max (eV)</th><th>Stokes shift (eV)</th><th>longest_bond_indices</th><th>Longest conjugated bond length</th></tr><tr><td>i64</td><td>str</td><td>str</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>str</td><td>f64</td><td>f64</td><td>f64</td><td>list[i64]</td><td>i64</td></tr></thead><tbody><tr><td>15095</td><td>&quot;c1ccc(C(=C(c2ccccc2)c2ccc(-c3c…</td><td>&quot;Cc1ccccc1&quot;</td><td>530.0</td><td>637.0</td><td>NaN</td><td>0.252</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1762.28448</td><td>&quot;https://doi.org/10.1016/j.dyep…</td><td>2.339282</td><td>1.946341</td><td>0.392941</td><td>[0, 1, … 157]</td><td>158</td></tr><tr><td>15096</td><td>&quot;c1ccc(C(=C(c2ccccc2)c2ccc(-c3c…</td><td>&quot;ClCCl&quot;</td><td>520.0</td><td>647.0</td><td>6.7</td><td>0.093</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1762.28448</td><td>&quot;https://doi.org/10.1016/j.dyep…</td><td>2.384268</td><td>1.916259</td><td>0.468009</td><td>[0, 1, … 157]</td><td>158</td></tr><tr><td>5658</td><td>&quot;c1ccc(-n2c(-c3ccc(N(c4ccc(-c5c…</td><td>&quot;ClCCl&quot;</td><td>376.0</td><td>432.0</td><td>1.68</td><td>0.206</td><td>NaN</td><td>NaN</td><td>NaN</td><td>60.4</td><td>62.7</td><td>1657.99344</td><td>&quot;https://doi.org/10.1016/j.orge…</td><td>3.297392</td><td>2.869952</td><td>0.42744</td><td>[0, 1, … 154]</td><td>155</td></tr><tr><td>5659</td><td>&quot;Cc1ccc(-n2c(-c3ccc(N(c4ccc(-c5…</td><td>&quot;ClCCl&quot;</td><td>376.0</td><td>435.0</td><td>1.68</td><td>0.282</td><td>NaN</td><td>NaN</td><td>NaN</td><td>61.2</td><td>62.1</td><td>1714.10096</td><td>&quot;https://doi.org/10.1016/j.orge…</td><td>3.297392</td><td>2.850159</td><td>0.447232</td><td>[1, 2, … 158]</td><td>155</td></tr><tr><td>5660</td><td>&quot;CC(C)(C)c1ccc(-n2c(-c3ccc(N(c4…</td><td>&quot;ClCCl&quot;</td><td>377.0</td><td>430.0</td><td>3.8</td><td>0.451</td><td>NaN</td><td>NaN</td><td>NaN</td><td>69.9</td><td>57.5</td><td>1882.42352</td><td>&quot;https://doi.org/10.1016/j.orge…</td><td>3.288645</td><td>2.883301</td><td>0.405345</td><td>[4, 5, … 170]</td><td>155</td></tr></tbody></table></div>



Checking the results, we find that the longest conjugated bond chain length in the optical dataset is 158 bonds! That's more than seven times longer than <i>beta-</i>carotene's. So it certainly makes sense that some of these molecules absorb visible light. For example, the molecule with the longest conjugated bond chain has its absorption maximum at 530 nm in the solvent of toluene.

## Checking for a correlation of absorption wavelength against longest conjugated bond chain length

Let's use [Polars' `plot`](https://docs.pola.rs/api/python/stable/reference/dataframe/plot.html) capability to plot absorption max against longest bond length to check for any trends. Altair can plot a maximum of 5,000 data points, and there are ~20,000 points in the optical dataset. We could enable [VegaFusion](https://vegafusion.io/) to allow for more points, but that would be an additional dependency and it seems not to work well with Google Colab. Instead, let's filter down to the ~7k unique chromophores in the dataset and plot the first 5k. Polars will select one solvent essentially at random for each chromophore, and the solvatochromic shift is generally not huge compared to the range we'll be plotting (from <200 to 850 nm), so that's a reasonable sampling. And 5k points should be plenty to discern if there's a trend.


```python
df_unique_chromophore = df.unique(subset="Chromophore")
```

We also need to drop the column `longest_bond_indices` because large lists of integers are not supported (we're not plotting them anyway)

Now we can make a plot of absorption maximum against longest conjugated bond length for the first 5k unique chromophores.


```python
df_unique_chromophore.slice(0, 5000).drop("longest_bond_indices").plot.scatter(
    x="Longest conjugated bond length", y="Absorption max (nm)"
)
```





<style>
  #altair-viz-e6c5a851cb5748afa3502008d54a0b84.vega-embed {
    width: 100%;
    display: flex;
  }

  #altair-viz-e6c5a851cb5748afa3502008d54a0b84.vega-embed details,
  #altair-viz-e6c5a851cb5748afa3502008d54a0b84.vega-embed details summary {
    position: relative;
  }
</style>
<div id="altair-viz-e6c5a851cb5748afa3502008d54a0b84"></div>
<script type="text/javascript">
  var VEGA_DEBUG = (typeof VEGA_DEBUG == "undefined") ? {} : VEGA_DEBUG;
  (function(spec, embedOpt){
    let outputDiv = document.currentScript.previousElementSibling;
    if (outputDiv.id !== "altair-viz-e6c5a851cb5748afa3502008d54a0b84") {
      outputDiv = document.getElementById("altair-viz-e6c5a851cb5748afa3502008d54a0b84");
    }
    const paths = {
      "vega": "https://cdn.jsdelivr.net/npm/vega@5?noext",
      "vega-lib": "https://cdn.jsdelivr.net/npm/vega-lib?noext",
      "vega-lite": "https://cdn.jsdelivr.net/npm/vega-lite@5.20.1?noext",
      "vega-embed": "https://cdn.jsdelivr.net/npm/vega-embed@6?noext",
    };

    function maybeLoadScript(lib, version) {
      var key = `${lib.replace("-", "")}_version`;
      return (VEGA_DEBUG[key] == version) ?
        Promise.resolve(paths[lib]) :
        new Promise(function(resolve, reject) {
          var s = document.createElement('script');
          document.getElementsByTagName("head")[0].appendChild(s);
          s.async = true;
          s.onload = () => {
            VEGA_DEBUG[key] = version;
            return resolve(paths[lib]);
          };
          s.onerror = () => reject(`Error loading script: ${paths[lib]}`);
          s.src = paths[lib];
        });
    }

    function showError(err) {
      outputDiv.innerHTML = `<div class="error" style="color:red;">${err}</div>`;
      throw err;
    }

    function displayChart(vegaEmbed) {
      vegaEmbed(outputDiv, spec, embedOpt)
        .catch(err => showError(`Javascript Error: ${err.message}<br>This usually means there's a typo in your chart specification. See the javascript console for the full traceback.`));
    }

    if(typeof define === "function" && define.amd) {
      requirejs.config({paths});
      require(["vega-embed"], displayChart, err => showError(`Error loading script: ${err.message}`));
    } else {
      maybeLoadScript("vega", "5")
        .then(() => maybeLoadScript("vega-lite", "5.20.1"))
        .then(() => maybeLoadScript("vega-embed", "6"))
        .catch(showError)
        .then(() => displayChart(vegaEmbed));
    }
</script>



There's not much of a correlation, presumably because of the varied molecular structures.

### Seeking a series of related molecules

Let's remove other variables by finding a series of molecules with a similar structure, but where the longest conjugated bond chain length increases. Given our dataset, let's start with molecules with a large number of connected conjugated bonds, then hopefully find a simple, linear molecule with a repeat unit that we can look for molecular analogues with fewer repeat units. We'll plot the top 50 molecules, which is the RDKit's default maximum for [`MolsToGridImage`](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.html#rdkit.Chem.Draw.MolsToGridImage).


```python
# Filter down to unique chromophores to avoid repeats in the molecular grid
df_unique_chromophore = df_unique_chromophore.sort(
    "Longest conjugated bond length", descending=True
)

# Extract columns as lists so we can plot molecules in a grid
unique_chromophores = df_unique_chromophore["Chromophore"].to_list()
tags = df_unique_chromophore["Tag"].to_list()
longest_conjugated_bond_lengths = df_unique_chromophore["Longest conjugated bond length"].to_list()
legends = [
    f"{longest_bond_length} bonds: tag {tag}"
    for (tag, longest_bond_length) in zip(tags, longest_conjugated_bond_lengths)
]
mols = [Chem.MolFromSmiles(chromophore) for chromophore in unique_chromophores]
Draw.MolsToGridImage(
    mols=mols,
    legends=legends,
    molsPerRow=5,
)
```




    
![Gallery of 50 molecules with the longest conjugated bond chain in the optical dataset](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_75_0.png)
    



While there are a lot of interesting-looking molecules, tag 19444 has a simple linear structure with a clear repeat unit that includes an anthracene moiety. This molecule has six of those repeat units, so we'll refer to it as the `n = 6` molecule. Let's show the molecule with its conjugated chain highlighted.


```python
size_wide_mol = (1000, 200)
df_n6 = df_unique_chromophore.filter(Tag=19444)
sml_n6 = df_n6[0]["Chromophore"].item()
mol_n6 = Chem.MolFromSmiles(sml_n6)
tags = df_n6[0]["Tag"].item()
longest_bond_length = df_n6[0]["Longest conjugated bond length"].item()
legend = f"{longest_bond_length} bonds: tag {tags}"
conjugated_bonds_n6 = df_n6[0]["longest_bond_indices"].item().to_list()
Draw.MolToImage(
    mol=mol_n6,
    size=size_wide_mol,
    legend=legend,
    highlightBonds=conjugated_bonds_n6,
)
```




    
![n=6 oligomer with conjugated bond chain highlighted](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_77_0.png)
    



By the way, sometimes it's difficult to see that all the specified bonds have been highlighted in the molecular image; there is in fact a contiguous chain from the first triple bond to the last.

Now let's show the repeat unit in the same orientation as the molecule, which requires rotating the atoms of the repeat unit by 90 degrees. We'll rotate serval molecules in this post by 90 degrees, so let's define a function to do that.


```python
def rotate_90(sml: str) -> Mol:
    """
    Rotates a molecule by 90 degrees based on its SMILES string and visualizes it.

    :param sml: SMILES string of the molecule to rotate.
    :returns: The rotated molecule.
    """

    mol = Chem.MolFromSmiles(sml)

    # Generate 2D coordinates if not already present
    AllChem.Compute2DCoords(mol)

    # Define the rotation matrix for 90 degrees (π/2 radians)
    theta = np.pi / 2  # 90 degrees in radians
    rotation_matrix = np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    )

    # Get the conformation of the molecule (the 2D coordinates)
    conf = mol.GetConformer()

    # Apply the rotation to each atom
    for i in range(mol.GetNumAtoms()):
        # Get the current x, y, z coordinates (z will stay 0.0 for 2D)
        pos = conf.GetAtomPosition(i)
        x, y = pos.x, pos.y  # Extract x, y coordinates

        # Apply the 90-degree rotation
        new_x, new_y = np.dot(rotation_matrix, [x, y])

        # Set the new coordinates, keeping z=0.0
        conf.SetAtomPosition(i, (new_x, new_y, 0.0))

    # Return the rotated molecule
    return mol
```

Now we can show the repeat unit in the same orientation that it appears in the `n = 6` molecule. The RDKit tends to recalculate molecular coordinates, so to be safe we assign the output image of `Draw.MolToImage` to a variable, then `display` that image.


```python
repeat_unit = "CC1=C2C=CC=CC2=C(C#C)C2=CC=CC=C12"
repeat_unit_mol = Chem.MolFromSmiles(repeat_unit)

repeat_unit_mol = rotate_90(repeat_unit)
image = Draw.MolToImage(
    repeat_unit_mol,
)
display(image)
```


    
![Anthracene-containing repeat unit. Off the left of the central ring are a carbon-carbon single bond and then a carbon-carbon triple bond; and off the righ of the central ring is carbon-carbon single bond.](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_82_0.png)
    


Let's check the dataframe for molecules with the same overall structure, with fewer repeat units. To do that, we'll need to check for a substructure in a molecule.

#### Filter to molecules in the series

Let's define a function to check how many of a given substructure there are in a molecule.


```python
def match_counts(
    sml: str,
    smls_to_match: Iterable[str] = None,
) -> Dict[str, int]:
    """
    Convert target and to-match SMILES to RDKit molecules, then count how many times the to-match molecules occur in the target molecules.

    :param sml: SMILES to convert to a molecule
    :param smls_to_match: One or more SMILES to find as a substructure
    :returns: Dictionary with name:count entry for each SMILES to match
    """
    mol = Chem.MolFromSmiles(sml)
    return_dict = dict()
    for name, sml_to_match in smls_to_match.items():
        mol_to_match = Chem.MolFromSmiles(sml_to_match)
        matches = mol.GetSubstructMatches(mol_to_match)
        return_dict[f"{name}_match_count"] = len(matches)
    return return_dict
```

To prevent double-counting of the number of repeat units, we use a substructure one atom longer than we showed above.


```python
smls_to_match = {
    # Repeat unit containing an anthracene moiety, with two triple bonds
    "repeat_unit": "C#CC1=C2C=CC=CC2=C(C#C)C2=CC=CC=C12",
}
```

Let's highlight one instance of the repeat unit, and of the terminal group, in the `n = 6` molecule. Because we'll search by atoms (the substructure), and then want to highlight the bonds as well, we'll define a function to get the bond indices connecting an iterable of atoms.


```python
def get_bond_indices_connecting_atoms(
    mol: Mol, atom_indices: Iterable[int]
) -> List[int]:
    """
    Given an RDKit molecule and a list of atom indices, return the bond indices
    that connect the atoms in the list.

    :param mol: RDKit molecule object
    :param atom_indices: List of atom indices to check for connectivity
    :return: List of bond indices that connect the atoms in atom_indices
    """
    bond_indices = []
    atom_set = set(atom_indices)  # Using a set for faster lookup

    # Iterate through all bonds in the molecule
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()

        # Check if both atoms of the bond are in the provided list of atom indices
        if begin_atom in atom_set and end_atom in atom_set:
            bond_indices.append(bond.GetIdx())

    return bond_indices
```


```python
# Get the first repeat unit substructure match
match_repeat_unit = mol_n6.GetSubstructMatches(repeat_unit_mol)[0]
bond_indices_repeat_unit = get_bond_indices_connecting_atoms(mol_n6, match_repeat_unit)
Draw.MolToImage(
    mol_n6,
    highlightAtoms=match_repeat_unit,
    highlightBonds=bond_indices_repeat_unit,
    size=size_wide_mol,
)
```




    
![n=6 oligomer where one anthracence-containing repeat unit is highlighted](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_91_0.png)
    



To avoid duplicate-column errors when the code is run more than once, we'll drop any existing match count columns. Then we'll count the occurrences of each substructure in each molecule.


```python
def add_match_counts(
    df: pl.DataFrame,
    smls_to_match: Iterable[str],
) -> pl.DataFrame:
    """
    For a dataframe with SMILES, add the count of matching substructures for SMILES to match.

    :param df: Polars dataframe of molecules with SMILES in Chromophore column
    :param smls_to_match: Iterable of SMILES for substructures
    :returns: Polars dataframe with match count columns added
    """
    df = df.drop(cs.ends_with("_match_count"))
    df = (
        df.with_columns(
            substructure_counts=pl.col("Chromophore").map_elements(
                lambda sml: match_counts(sml, smls_to_match), return_dtype=pl.Struct
            )
        )
        .unnest("substructure_counts")
        .sort("Longest conjugated bond length", descending=True)
        .shrink_to_fit()
    )
    return df
```


```python
df = add_match_counts(
    df=df,
    smls_to_match=smls_to_match,
)
```

Let's filter down to molecules that contain the repeat unit by requiring the repeat unit match count be greater than zero.


```python
df.filter((pl.col("repeat_unit_match_count") > 0)).sort(
    pl.col("repeat_unit_match_count")
)
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (51, 20)</small><table border="1" class="dataframe"><thead><tr><th>Tag</th><th>Chromophore</th><th>Solvent</th><th>Absorption max (nm)</th><th>Emission max (nm)</th><th>Lifetime (ns)</th><th>Quantum yield</th><th>log(e/mol-1 dm3 cm-1)</th><th>abs FWHM (cm-1)</th><th>emi FWHM (cm-1)</th><th>abs FWHM (nm)</th><th>emi FWHM (nm)</th><th>Molecular weight (g mol-1)</th><th>Reference</th><th>Absorption max (eV)</th><th>Emission max (eV)</th><th>Stokes shift (eV)</th><th>longest_bond_indices</th><th>Longest conjugated bond length</th><th>repeat_unit_match_count</th></tr><tr><td>i64</td><td>str</td><td>str</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>str</td><td>f64</td><td>f64</td><td>f64</td><td>list[i64]</td><td>i64</td><td>i64</td></tr></thead><tbody><tr><td>18923</td><td>&quot;CCCCCCCCn1c(=O)c2cc3c(C#Cc4ccc…</td><td>&quot;ClC(Cl)Cl&quot;</td><td>582.0</td><td>641.0</td><td>3.4</td><td>0.53</td><td>4.518514</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1075.36204</td><td>&quot;DOI: 10.1021/acs.joc.6b00364&quot;</td><td>2.130274</td><td>1.934196</td><td>0.196078</td><td>[8, 9, … 91]</td><td>76</td><td>1</td></tr><tr><td>19570</td><td>&quot;O=c1c2ccccc2c2nc3cc4c(C#C[Si](…</td><td>&quot;ClCCl&quot;</td><td>731.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1239.80484</td><td>&quot;DOI: 10.1021/acs.joc.9b02756&quot;</td><td>1.696059</td><td>NaN</td><td>NaN</td><td>[0, 1, … 107]</td><td>64</td><td>1</td></tr><tr><td>19577</td><td>&quot;O=c1c2ccccc2c2nc3cc4c(C#C[Si](…</td><td>&quot;O=c1c2ccccc2c2nc3cc4c(C#C[Si](…</td><td>743.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1239.80484</td><td>&quot;DOI: 10.1021/acs.joc.9b02756&quot;</td><td>1.668667</td><td>NaN</td><td>NaN</td><td>[0, 1, … 107]</td><td>64</td><td>1</td></tr><tr><td>19367</td><td>&quot;COc1ccc(C#Cc2c3cc4ccccc4cc3c(C…</td><td>&quot;ClCCl&quot;</td><td>666.0</td><td>706.0</td><td>4.6</td><td>0.03</td><td>4.436163</td><td>NaN</td><td>NaN</td><td>NaN</td><td>62.3</td><td>538.64444</td><td>&quot;DOI: 10.1021/jo501696d&quot;</td><td>1.861591</td><td>1.756118</td><td>0.105473</td><td>[1, 2, … 47]</td><td>46</td><td>1</td></tr><tr><td>19516</td><td>&quot;COc1c2ccccc2c(OC)c2c(C#Cc3cccc…</td><td>&quot;ClCCl&quot;</td><td>704.0</td><td>718.0</td><td>NaN</td><td>0.11</td><td>4.130334</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>538.64444</td><td>&quot;DOI: 10.1021/jo0710331&quot;</td><td>1.761107</td><td>1.726768</td><td>0.034339</td><td>[1, 2, … 47]</td><td>46</td><td>1</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>19440</td><td>&quot;CC(C)[Si](C#Cc1c2ccc(C(C)(C)C)…</td><td>&quot;Cc1ccccc1&quot;</td><td>523.0</td><td>541.0</td><td>NaN</td><td>0.02</td><td>4.694605</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>963.6346</td><td>&quot;DOI: 10.1021/acs.joc.8b00311&quot;</td><td>2.370591</td><td>2.291718</td><td>0.078874</td><td>[4, 5, … 74]</td><td>39</td><td>2</td></tr><tr><td>19441</td><td>&quot;CC(C)[Si](C#Cc1c2ccc(C(C)(C)C)…</td><td>&quot;Cc1ccccc1&quot;</td><td>570.0</td><td>589.0</td><td>NaN</td><td>0.022</td><td>4.781037</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1276.08916</td><td>&quot;DOI: 10.1021/acs.joc.8b00311&quot;</td><td>2.175122</td><td>2.104956</td><td>0.070165</td><td>[4, 5, … 101]</td><td>58</td><td>3</td></tr><tr><td>19442</td><td>&quot;CC(C)[Si](C#Cc1c2ccccc2c(C#Cc2…</td><td>&quot;Cc1ccccc1&quot;</td><td>582.0</td><td>609.0</td><td>NaN</td><td>0.02</td><td>4.975432</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1139.68356</td><td>&quot;DOI: 10.1021/acs.joc.8b00311&quot;</td><td>2.130274</td><td>2.035828</td><td>0.094446</td><td>[4, 5, … 96]</td><td>77</td><td>4</td></tr><tr><td>19443</td><td>&quot;CC(C)[Si](C#Cc1c2ccccc2c(C#Cc2…</td><td>&quot;Cc1ccccc1&quot;</td><td>585.0</td><td>623.0</td><td>NaN</td><td>0.018</td><td>5.130334</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1339.92308</td><td>&quot;DOI: 10.1021/acs.joc.8b00311&quot;</td><td>2.119349</td><td>1.990079</td><td>0.12927</td><td>[4, 5, … 115]</td><td>96</td><td>5</td></tr><tr><td>19444</td><td>&quot;CC(C)[Si](C#Cc1c2ccccc2c(C#Cc2…</td><td>&quot;Cc1ccccc1&quot;</td><td>589.0</td><td>629.0</td><td>NaN</td><td>NaN</td><td>5.178977</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1540.1626</td><td>&quot;DOI: 10.1021/acs.joc.8b00311&quot;</td><td>2.104956</td><td>1.971096</td><td>0.133861</td><td>[4, 5, … 134]</td><td>115</td><td>6</td></tr></tbody></table></div>



That's too many molecules--we expect as many as six, for `n = 1-6`. Let's examine the molecular structures to understand why.

Because we'll plot several instances of molecular grids from subsets of the dataframe, we define a function to plot first ten molecules and legends based on the columns in the dataframe.


```python
def df_to_grid_image(
    df: pl.DataFrame,
    subImgSize: Tuple[int, int] = (200, 200),
    molsPerRow=3,
) -> Image.Image:
    """
    Create a molecular grid image from a Polars dataframe.

    :param df: Polars dataframe
    :param subImgSize: The size for each molecule
    :param molsPerRow: The number of molecules in each row
    :returns: Molecular grid image
    """
    matching = df["Chromophore"].to_list()
    mols = [Chem.MolFromSmiles(match) for match in matching]
    nums_match = df["repeat_unit_match_count"].to_list()
    absorbances_nm = df["Absorption max (nm)"].to_list()
    tags = df["Tag"].to_list()
    legends = []
    for index, num_match in enumerate(nums_match):
        this_tag = tags[index]
        legend = f"Tag {this_tag}: {num_match} unit(s)"
        abs_nm = absorbances_nm[index]
        if not math.isnan(abs_nm):
            legend += f" {abs_nm}nm"
        legends.append(legend)
    conjugated_bonds = [get_longest_conjugated_bond_chain(mol) for mol in mols]
    print(f"{len(mols)} molecules")

    dwg = Draw.MolsToGridImage(
        mols=mols,
        legends=legends,
        maxMols=10,
        molsPerRow=molsPerRow,
        highlightBondLists=conjugated_bonds,
        subImgSize=subImgSize,
    )
    return dwg
```


```python
df_to_grid_image(
    df.filter((pl.col("repeat_unit_match_count") > 0)).sort(
        pl.col("repeat_unit_match_count")
    )
)
```

    51 molecules





    
![Grid of 10 molecules with the anthracene-containing repeat unit](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_100_1.png)
    



None of the first 10 molecules have the two terminal triisopropylsyl groups. Let's specify that we want two such terminal groups by adding the terminal group to the list of SMILES to match, then filtering to molecules with two of those terminal groups.


```python
smls_to_match = smls_to_match | {
    # Terminal triisopropylsyl group
    "terminal": "CC(C)[SiH](C(C)C)C(C)C"
}
```

Here are the substructures we're searching for now:


```python
mols_to_match = []
for sml in smls_to_match.values():
    mol = rotate_90(sml)
    mols_to_match.append(mol)
image = Draw.MolsToGridImage(
    mols=mols_to_match,
    legends=smls_to_match.keys(),
    molsPerRow=4,
    subImgSize=(300, 350),
)
display(image)
```


    
![Anthracene-containing repeat unit, and triisopropylsyl group](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_104_0.png)
    



```python
df = add_match_counts(
    df=df,
    smls_to_match=smls_to_match,
)
```


```python
df_to_grid_image(
    df.filter(
        (pl.col("repeat_unit_match_count") > 0) & (pl.col("terminal_match_count") == 2)
    ).sort(pl.col("repeat_unit_match_count"))
)
```

    19 molecules





    
![Grid of 10 molecules with the anthracene-containing repeat unit and two triisopropylsyl groups](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_106_1.png)
    



That helps, but we're getting a lot of fused-ring systems with more than three rings. Let's exclude those by adding SMILES to match that contain four fused rings, where the fourth ring is either all carbons, or contains two nitrogens. We'll then filter to molecules whose match count is zero four-fused-ring substructures.


```python
smls_to_match = smls_to_match | {
    # Four benzene rings
    "too_many_rings": "C#CC1=C2C=C3C=CC=CC3=CC2=C(C#C)C2=CC=CC=C12",
    # Four rings: Three benzene and one dinitrogen ring
    "too_many_rings_N": "C#CC1=C2C=C3N=CC=NC3=CC2=C(C#C)C2=CC=CC=C12",
}
```

Here are all the substructures we're searching for now, where the first row are substructures we want to include, while the second row is substructures we want to exclude.


```python
mols_to_match = []
for sml in smls_to_match.values():
    mol = rotate_90(sml)
    mols_to_match.append(mol)
image = Draw.MolsToGridImage(
    mols=mols_to_match,
    legends=smls_to_match.keys(),
    molsPerRow=2,
    subImgSize=(300, 350),
)
display(image)
```


    
![Top row: Anthracene-containing repeat unit, and triisopropylsyl group. Bottom row: Similar repeat unit but with four fused rings instead of three, with the fourth ring either all carbons or with two nitrogens.](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_110_0.png)
    



```python
df = add_match_counts(
    df=df,
    smls_to_match=smls_to_match,
)
```

Next we exclude the four-fused-ring substructures by setting their count to zero.


```python
df_matches = df.filter(
    (pl.col("repeat_unit_match_count") > 0)
    & (pl.col("terminal_match_count") == 2)
    & (pl.col("too_many_rings_match_count") == 0)
    & (pl.col("too_many_rings_N_match_count") == 0)
).sort(pl.col("repeat_unit_match_count"))
df_to_grid_image(
    df_matches,
    subImgSize=(900, 300),
    molsPerRow=3,
)
```

    6 molecules





    
![n=1-6 oligomers](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_113_1.png)
    




```python
df_matches
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (6, 23)</small><table border="1" class="dataframe"><thead><tr><th>Tag</th><th>Chromophore</th><th>Solvent</th><th>Absorption max (nm)</th><th>Emission max (nm)</th><th>Lifetime (ns)</th><th>Quantum yield</th><th>log(e/mol-1 dm3 cm-1)</th><th>abs FWHM (cm-1)</th><th>emi FWHM (cm-1)</th><th>abs FWHM (nm)</th><th>emi FWHM (nm)</th><th>Molecular weight (g mol-1)</th><th>Reference</th><th>Absorption max (eV)</th><th>Emission max (eV)</th><th>Stokes shift (eV)</th><th>longest_bond_indices</th><th>Longest conjugated bond length</th><th>repeat_unit_match_count</th><th>terminal_match_count</th><th>too_many_rings_match_count</th><th>too_many_rings_N_match_count</th></tr><tr><td>i64</td><td>str</td><td>str</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>f64</td><td>str</td><td>f64</td><td>f64</td><td>f64</td><td>list[i64]</td><td>i64</td><td>i64</td><td>i64</td><td>i64</td><td>i64</td></tr></thead><tbody><tr><td>19439</td><td>&quot;CC(C)[Si](C#Cc1c2ccc(C(C)(C)C)…</td><td>&quot;Cc1ccccc1&quot;</td><td>440.0</td><td>461.0</td><td>NaN</td><td>0.92</td><td>4.4133</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>651.18004</td><td>&quot;DOI: 10.1021/acs.joc.8b00311&quot;</td><td>2.817771</td><td>2.689413</td><td>0.128358</td><td>[4, 5, … 47]</td><td>20</td><td>1</td><td>2</td><td>0</td><td>0</td></tr><tr><td>19440</td><td>&quot;CC(C)[Si](C#Cc1c2ccc(C(C)(C)C)…</td><td>&quot;Cc1ccccc1&quot;</td><td>523.0</td><td>541.0</td><td>NaN</td><td>0.02</td><td>4.694605</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>963.6346</td><td>&quot;DOI: 10.1021/acs.joc.8b00311&quot;</td><td>2.370591</td><td>2.291718</td><td>0.078874</td><td>[4, 5, … 74]</td><td>39</td><td>2</td><td>2</td><td>0</td><td>0</td></tr><tr><td>19441</td><td>&quot;CC(C)[Si](C#Cc1c2ccc(C(C)(C)C)…</td><td>&quot;Cc1ccccc1&quot;</td><td>570.0</td><td>589.0</td><td>NaN</td><td>0.022</td><td>4.781037</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1276.08916</td><td>&quot;DOI: 10.1021/acs.joc.8b00311&quot;</td><td>2.175122</td><td>2.104956</td><td>0.070165</td><td>[4, 5, … 101]</td><td>58</td><td>3</td><td>2</td><td>0</td><td>0</td></tr><tr><td>19442</td><td>&quot;CC(C)[Si](C#Cc1c2ccccc2c(C#Cc2…</td><td>&quot;Cc1ccccc1&quot;</td><td>582.0</td><td>609.0</td><td>NaN</td><td>0.02</td><td>4.975432</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1139.68356</td><td>&quot;DOI: 10.1021/acs.joc.8b00311&quot;</td><td>2.130274</td><td>2.035828</td><td>0.094446</td><td>[4, 5, … 96]</td><td>77</td><td>4</td><td>2</td><td>0</td><td>0</td></tr><tr><td>19443</td><td>&quot;CC(C)[Si](C#Cc1c2ccccc2c(C#Cc2…</td><td>&quot;Cc1ccccc1&quot;</td><td>585.0</td><td>623.0</td><td>NaN</td><td>0.018</td><td>5.130334</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1339.92308</td><td>&quot;DOI: 10.1021/acs.joc.8b00311&quot;</td><td>2.119349</td><td>1.990079</td><td>0.12927</td><td>[4, 5, … 115]</td><td>96</td><td>5</td><td>2</td><td>0</td><td>0</td></tr><tr><td>19444</td><td>&quot;CC(C)[Si](C#Cc1c2ccccc2c(C#Cc2…</td><td>&quot;Cc1ccccc1&quot;</td><td>589.0</td><td>629.0</td><td>NaN</td><td>NaN</td><td>5.178977</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>1540.1626</td><td>&quot;DOI: 10.1021/acs.joc.8b00311&quot;</td><td>2.104956</td><td>1.971096</td><td>0.133861</td><td>[4, 5, … 134]</td><td>115</td><td>6</td><td>2</td><td>0</td><td>0</td></tr></tbody></table></div>



Great! We have just the molecules we want, and the dataset contains all six oligomers with the anthracene repeat unit: `n = 1-6`.

By the way, those four substructure criteria are sufficient to get just the set of molecules we want from this dataset. If we had a larger dataset, we might have to add more criteria to exclude other molecules.

Also, the solvent is toluene in all cases, meaning there are no different solvent effects to consider.


```python
Chem.MolFromSmiles(df_matches[0]["Solvent"].item())
```




    
![Toleuene molecular structure](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_118_0.png)
    



It turns out these molecules all come from one paper, "Synthesis and Electronic Properties of Length-Defined 9,10-Anthrylene−Butadiynylene Oligomers" by Nagaoka et al., DOI [10.1021/acs.joc.8b00311](https://pubs.acs.org/doi/full/10.1021/acs.joc.8b00311). They actually made this series of molecules for their optical properties: "These oligomers will be attractive as core units in the molecular design of linear π-conjugated compounds for functional dyes and electronic devices." Their calculated Stokes shift are the same as in our table above.

They also note that the "reaction mixture instantly turned deep purple as a result of the formation of a mixture of oligomers." This is in agreement with our earlier statement about complimentary colors: Because the `n = 1-6` oligomers absorb from 440-589 nm, the complimentary colors will include purple.

There appear to be some minor discrepancies between the molecular structures in the database and Nagaoka's paper: all of the anthracene units should have two <i>tert-</i>butyl groups, and their orientation should alternate. We manually updated the SMILES below. These changes don't affect which molecules were filtered by our four criteria.


```python
sml_n6_corrected = "CC(C)[Si](C#Cc%17c1cc(C(C)(C)C)ccc1c(C#Cc%15c2ccc(C(C)(C)C)cc2c(C#Cc%13c3cc(C(C)(C)C)ccc3c(C#Cc%11c4ccc(C(C)(C)C)cc4c(C#Cc9c5cc(C(C)(C)C)ccc5c(C#Cc7c6ccc(C(C)(C)C)cc6c(C#C[Si](C(C)C)(C(C)C)C(C)C)c8ccc(C(C)(C)C)cc78)c%10cc(C(C)(C)C)ccc9%10)c%12ccc(C(C)(C)C)cc%11%12)c%14cc(C(C)(C)C)ccc%13%14)c%16ccc(C(C)(C)C)cc%15%16)c%18cc(C(C)(C)C)ccc%17%18)(C(C)C)C(C)C"
sml_n5_corrected = "CC(C)[Si](C#Cc%15c2ccc(C(C)(C)C)cc2c(C#Cc%13c3cc(C(C)(C)C)ccc3c(C#Cc%11c4ccc(C(C)(C)C)cc4c(C#Cc9c5cc(C(C)(C)C)ccc5c(C#Cc7c6ccc(C(C)(C)C)cc6c(C#C[Si](C(C)C)(C(C)C)C(C)C)c8ccc(C(C)(C)C)cc78)c%10cc(C(C)(C)C)ccc9%10)c%12ccc(C(C)(C)C)cc%11%12)c%14cc(C(C)(C)C)ccc%13%14)c%16ccc(C(C)(C)C)cc%15%16)(C(C)C)C(C)C"
sml_n4_corrected = "CC(C)[Si](C#Cc%11c1cc(C(C)(C)C)ccc1c(C#Cc9c2ccc(C(C)(C)C)cc2c(C#Cc7c3cc(C(C)(C)C)ccc3c(C#Cc5c4ccc(C(C)(C)C)cc4c(C#C[Si](C(C)C)(C(C)C)C(C)C)c6ccc(C(C)(C)C)cc56)c8cc(C(C)(C)C)ccc78)c%10ccc(C(C)(C)C)cc9%10)c%12cc(C(C)(C)C)ccc%11%12)(C(C)C)C(C)C"
sml_n3_corrected = "CC(C)[Si](C#Cc8c1cc(C(C)(C)C)ccc1c(C#Cc6c2ccc(C(C)(C)C)cc2c(C#Cc4c3cc(C(C)(C)C)ccc3c(C#C[Si](C(C)C)(C(C)C)C(C)C)c5cc(C(C)(C)C)ccc45)c7ccc(C(C)(C)C)cc67)c9cc(C(C)(C)C)ccc89)(C(C)C)C(C)C"
sml_n2_corrected = "CC(C)[Si](C#Cc5c1cc(C(C)(C)C)ccc1c(C#Cc3c2ccc(C(C)(C)C)cc2c(C#C[Si](C(C)C)(C(C)C)C(C)C)c4ccc(C(C)(C)C)cc34)c6cc(C(C)(C)C)ccc56)(C(C)C)C(C)C"

sml_n5 = df_unique_chromophore.filter(Tag=19443)[0]["Chromophore"].item()
sml_n4 = df_unique_chromophore.filter(Tag=19442)[0]["Chromophore"].item()
sml_n3 = df_unique_chromophore.filter(Tag=19441)[0]["Chromophore"].item()
sml_n2 = df_unique_chromophore.filter(Tag=19440)[0]["Chromophore"].item()

df_matches = (
    df_matches.with_columns(pl.col("Chromophore").replace(sml_n6, sml_n6_corrected))
    .with_columns(pl.col("Chromophore").replace(sml_n5, sml_n5_corrected))
    .with_columns(pl.col("Chromophore").replace(sml_n4, sml_n4_corrected))
    .with_columns(pl.col("Chromophore").replace(sml_n3, sml_n3_corrected))
    .with_columns(pl.col("Chromophore").replace(sml_n2, sml_n2_corrected))
)
```


```python
df_to_grid_image(
    df_matches,
    subImgSize=(900, 300),
    molsPerRow=3,
)
```

    6 molecules





    
![n=1-6 oligomers with corrected structures](/images/2024-10-15-Color-from-Conjugation_files/2024-10-15-Color-from-Conjugation_122_1.png)
    



As a side note, from examining the dataframe, the [photoluminescence](https://en.wikipedia.org/wiki/Photoluminescence) quantum yield is 0.92 for `n = 1`, and about 0.02 for `n = 2-5` (a value is not provided for `n = 6`). Nagaoka et al. attribute this to the fact that "the presence of diacetylene moieties facilitates quenching via nonradiative pathways, as reported for other butadiyne compounds."

## Plot results to check for correlation between color and conjugation chain length

Now that we have a consistent series of molecules, let's plot their optical properties to check for trends. We'll use the [altair](https://altair-viz.github.io/) library that provides Polars' plotting capability. 

We start by setting our shape scheme for the absorption, emission, and Stokes shift series.


```python
emission_shape = "circle"
absorption_shape = "triangle"
stokes_shape = "square"
```

Next we calculate plot ranges based on the data. For all axes, we'll allow for a slight buffer so our data points aren't on an edge of the plot.


```python
# Calculate the minimum and maximum y-value in nm across the absorption and emission series
min_y_value_nm = min(
    df_matches["Absorption max (nm)"].min(), df_matches["Emission max (nm)"].min()
)
y_min_nm = min_y_value_nm * 0.95
max_y_value_nm = max(
    df_matches["Absorption max (nm)"].max(), df_matches["Emission max (nm)"].max()
)
y_max_nm = max_y_value_nm * 1.05

# Calculate the minimum and maximum y-value in eV across the absorption and emission series
min_y_value_eV = min(
    df_matches["Absorption max (eV)"].min(), df_matches["Emission max (eV)"].min()
)
y_min_eV = min_y_value_eV * 0.95
max_y_value_eV = max(
    df_matches["Absorption max (eV)"].max(), df_matches["Emission max (eV)"].max()
)
y_max_eV = max_y_value_eV * 1.05

# Calculate the maximum y-value in eV for the Stokes shifts
stokes_y_max = df_matches["Stokes shift (eV)"].max() * 1.1

# Calculate the maximum x-value for all series, based on the maximum repeat unit count
max_repeat_count = df_matches["repeat_unit_match_count"].max() + 0.5
```

Now we can plot the data. Let's start with the wavelength units in the source dataset, namely nanometers (nm). We'll make a scatter plot for the absorption series, then one for the emission series, and layer them on top of each other. To emphasize the color of the light, let's make each symbol the approximate color of that light, for example 440 nm is violet. For the legend, we'll make a separate "plot" and concatenate it horizontally with the main plot.


```python
# Order of curves from top to bottom
shapes = [emission_shape, absorption_shape]

# Larger data point size for better visibility
point_size = 150

# Thickness of black outline on data points
stroke_width = 2

# Define a color scale that corresponds to the visible spectrum
color_scale = alt.Scale(
    domain=[380, 700],  # Wavelength range (approximate visible light)
    range=[
        "#8B00FF",  # Violet (~380 nm)
        "#0000FF",  # Blue (~450 nm)
        "#00FF00",  # Green (~520 nm)
        "#FFFF00",  # Yellow (~580 nm)
        "#FF7F00",  # Orange (~600 nm)
        "#FF0000",  # Red (~700 nm)
    ],
)

# Absorption scatter plot with black stroke and fill color based on wavelength
scatter_absorption = (
    alt.Chart(df_matches)
    # Black border with a fill color
    .mark_point(size=point_size, filled=True, stroke="black", strokeWidth=stroke_width)
    .encode(
        x="repeat_unit_match_count",
        y="Absorption max (nm)",
        color=alt.Color(
            "Absorption max (nm)",
            # Fill color corresponds to wavelength
            scale=color_scale,
            title="Wavelength (nm)",
            legend=None,
        ),
        shape=alt.value(absorption_shape),
        tooltip=["repeat_unit_match_count", "Absorption max (nm)"],
    )
    .properties(title="Absorption max")
)

# Emission scatter plot with black stroke and fill color based on wavelength
scatter_emission = (
    alt.Chart(df_matches)
    # Black border with a fill color
    .mark_point(size=point_size, filled=True, stroke="black", strokeWidth=stroke_width)
    .encode(
        x="repeat_unit_match_count",
        y="Emission max (nm)",
        color=alt.Color(
            "Emission max (nm)",
            # Fill color corresponds to wavelength
            scale=color_scale,
            title="Wavelength (nm)",
            legend=None,
        ),
        # Shape for emission
        shape=alt.value(emission_shape),
        tooltip=["repeat_unit_match_count", "Emission max (nm)"],
    )
    .properties(title="Emission max")
)

# Manually create a legend for shapes with the reordered entries (Emission on top, Absorption on bottom)
legend_data = pl.DataFrame(
    {
        "Type": [
            "Emission max (eV)",
            "Absorption max (eV)",
        ],
        # Reordered: Emission first, Absorption second
        "Shape": shapes,
    }
)

legend = (
    alt.Chart(legend_data)
    .mark_point(size=point_size)
    .encode(
        y=alt.Y(
            "Type:N",
            sort=[
                "Emission max (eV)",
                "Absorption max (eV)",
            ],  # Explicit sort order
            axis=alt.Axis(orient="right"),
            title="",
        ),
        shape=alt.Shape(
            "Shape:N",
            scale=alt.Scale(domain=shapes, range=shapes),
            # Disable shape legend
            legend=None,
        ),
        # Keep legend shapes black for distinction
        color=alt.value("black"),
    )
)

# Combine the scatter plots
scatter_combined = (
    alt.layer(scatter_absorption, scatter_emission)
    .properties(title="Optical Properties of Oligomer Series")
    .encode(
        y=alt.Y(
            "Absorption max (nm)",
            scale=alt.Scale(domain=[y_min_nm, y_max_nm]),
            title="Absorption or Emission Max (nm)",
        ),
        x=alt.X(
            "repeat_unit_match_count",
            scale=alt.Scale(domain=[0.5, max_repeat_count]),
            title="Number of Repeat Units",
        ),
    )
    .interactive()
)

# Combine the scatter plot with the legend
final_chart = alt.hconcat(scatter_combined, legend)

# Display the chart
final_chart.show()
```



<style>
  #altair-viz-401332e1269c4cdca82c9ec3e64b9c60.vega-embed {
    width: 100%;
    display: flex;
  }

  #altair-viz-401332e1269c4cdca82c9ec3e64b9c60.vega-embed details,
  #altair-viz-401332e1269c4cdca82c9ec3e64b9c60.vega-embed details summary {
    position: relative;
  }
</style>
<div id="altair-viz-401332e1269c4cdca82c9ec3e64b9c60"></div>
<script type="text/javascript">
  var VEGA_DEBUG = (typeof VEGA_DEBUG == "undefined") ? {} : VEGA_DEBUG;
  (function(spec, embedOpt){
    let outputDiv = document.currentScript.previousElementSibling;
    if (outputDiv.id !== "altair-viz-401332e1269c4cdca82c9ec3e64b9c60") {
      outputDiv = document.getElementById("altair-viz-401332e1269c4cdca82c9ec3e64b9c60");
    }
    const paths = {
      "vega": "https://cdn.jsdelivr.net/npm/vega@5?noext",
      "vega-lib": "https://cdn.jsdelivr.net/npm/vega-lib?noext",
      "vega-lite": "https://cdn.jsdelivr.net/npm/vega-lite@5.20.1?noext",
      "vega-embed": "https://cdn.jsdelivr.net/npm/vega-embed@6?noext",
    };

    function maybeLoadScript(lib, version) {
      var key = `${lib.replace("-", "")}_version`;
      return (VEGA_DEBUG[key] == version) ?
        Promise.resolve(paths[lib]) :
        new Promise(function(resolve, reject) {
          var s = document.createElement('script');
          document.getElementsByTagName("head")[0].appendChild(s);
          s.async = true;
          s.onload = () => {
            VEGA_DEBUG[key] = version;
            return resolve(paths[lib]);
          };
          s.onerror = () => reject(`Error loading script: ${paths[lib]}`);
          s.src = paths[lib];
        });
    }

    function showError(err) {
      outputDiv.innerHTML = `<div class="error" style="color:red;">${err}</div>`;
      throw err;
    }

    function displayChart(vegaEmbed) {
      vegaEmbed(outputDiv, spec, embedOpt)
        .catch(err => showError(`Javascript Error: ${err.message}<br>This usually means there's a typo in your chart specification. See the javascript console for the full traceback.`));
    }

    if(typeof define === "function" && define.amd) {
      requirejs.config({paths});
      require(["vega-embed"], displayChart, err => showError(`Error loading script: ${err.message}`));
    } else {
      maybeLoadScript("vega", "5")
        .then(() => maybeLoadScript("vega-lite", "5.20.1"))
        .then(() => maybeLoadScript("vega-embed", "6"))
        .catch(showError)
        .then(() => displayChart(vegaEmbed));
    }
  })({"config": {"view": {"continuousWidth": 300, "continuousHeight": 300}}, "hconcat": [{"layer": [{"mark": {"type": "point", "filled": true, "size": 150, "stroke": "black", "strokeWidth": 2}, "encoding": {"color": {"field": "Absorption max (nm)", "legend": null, "scale": {"domain": [380, 700], "range": ["#8B00FF", "#0000FF", "#00FF00", "#FFFF00", "#FF7F00", "#FF0000"]}, "title": "Wavelength (nm)", "type": "quantitative"}, "shape": {"value": "triangle"}, "tooltip": [{"field": "repeat_unit_match_count", "type": "quantitative"}, {"field": "Absorption max (nm)", "type": "quantitative"}], "x": {"field": "repeat_unit_match_count", "type": "quantitative"}, "y": {"field": "Absorption max (nm)", "type": "quantitative"}}, "name": "view_1", "title": "Absorption max"}, {"mark": {"type": "point", "filled": true, "size": 150, "stroke": "black", "strokeWidth": 2}, "encoding": {"color": {"field": "Emission max (nm)", "legend": null, "scale": {"domain": [380, 700], "range": ["#8B00FF", "#0000FF", "#00FF00", "#FFFF00", "#FF7F00", "#FF0000"]}, "title": "Wavelength (nm)", "type": "quantitative"}, "shape": {"value": "circle"}, "tooltip": [{"field": "repeat_unit_match_count", "type": "quantitative"}, {"field": "Emission max (nm)", "type": "quantitative"}], "x": {"field": "repeat_unit_match_count", "type": "quantitative"}, "y": {"field": "Emission max (nm)", "type": "quantitative"}}, "title": "Emission max"}], "data": {"name": "data-fcb24888a4ab58ce1808fb0997ef05d4"}, "encoding": {"x": {"field": "repeat_unit_match_count", "scale": {"domain": [0.5, 6.5]}, "title": "Number of Repeat Units", "type": "quantitative"}, "y": {"field": "Absorption max (nm)", "scale": {"domain": [418.0, 660.45]}, "title": "Absorption or Emission Max (nm)", "type": "quantitative"}}, "title": "Optical Properties of Oligomer Series"}, {"data": {"name": "data-23d98a78cf3a4d67e6e8e922ac3376db"}, "mark": {"type": "point", "size": 150}, "encoding": {"color": {"value": "black"}, "shape": {"field": "Shape", "legend": null, "scale": {"domain": ["circle", "triangle"], "range": ["circle", "triangle"]}, "type": "nominal"}, "y": {"axis": {"orient": "right"}, "field": "Type", "sort": ["Emission max (eV)", "Absorption max (eV)"], "title": "", "type": "nominal"}}}], "params": [{"name": "param_2", "select": {"type": "interval", "encodings": ["x", "y"]}, "bind": "scales", "views": ["view_1"]}], "$schema": "https://vega.github.io/schema/vega-lite/v5.20.1.json", "datasets": {"data-fcb24888a4ab58ce1808fb0997ef05d4": [{"Tag": 19439, "Chromophore": "CC(C)[Si](C#Cc1c2ccc(C(C)(C)C)cc2c(C#C[Si](C(C)C)(C(C)C)C(C)C)c2ccc(C(C)(C)C)cc12)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 440.0, "Emission max (nm)": 461.0, "Lifetime (ns)": NaN, "Quantum yield": 0.92, "log(e/mol-1 dm3 cm-1)": 4.413299764, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 651.18004, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.817771188195776, "Emission max (eV)": 2.6894128477356647, "Stokes shift (eV)": 0.1283583404601112, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 29, 30, 31, 32, 37, 38, 45, 46, 47], "Longest conjugated bond length": 20, "repeat_unit_match_count": 1, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}, {"Tag": 19440, "Chromophore": "CC(C)[Si](C#Cc5c1cc(C(C)(C)C)ccc1c(C#Cc3c2ccc(C(C)(C)C)cc2c(C#C[Si](C(C)C)(C(C)C)C(C)C)c4ccc(C(C)(C)C)cc34)c6cc(C(C)(C)C)ccc56)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 523.0, "Emission max (nm)": 541.0, "Lifetime (ns)": NaN, "Quantum yield": 0.02, "log(e/mol-1 dm3 cm-1)": 4.694605199, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 963.6346, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.3705914393998877, "Emission max (eV)": 2.29171778707235, "Stokes shift (eV)": 0.07887365232753751, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 28, 29, 30, 31, 32, 43, 44, 45, 46, 51, 52, 53, 54, 55, 56, 61, 62, 69, 70, 71, 72, 73, 74], "Longest conjugated bond length": 39, "repeat_unit_match_count": 2, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}, {"Tag": 19441, "Chromophore": "CC(C)[Si](C#Cc8c1cc(C(C)(C)C)ccc1c(C#Cc6c2ccc(C(C)(C)C)cc2c(C#Cc4c3cc(C(C)(C)C)ccc3c(C#C[Si](C(C)C)(C(C)C)C(C)C)c5cc(C(C)(C)C)ccc45)c7ccc(C(C)(C)C)cc67)c9cc(C(C)(C)C)ccc89)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 570.0, "Emission max (nm)": 589.0, "Lifetime (ns)": NaN, "Quantum yield": 0.022, "log(e/mol-1 dm3 cm-1)": 4.781036939, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 1276.08916, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.175121618958143, "Emission max (eV)": 2.104956405443364, "Stokes shift (eV)": 0.070165213514779, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 42, 43, 44, 45, 46, 57, 58, 59, 60, 65, 66, 67, 68, 69, 70, 75, 76, 77, 78, 79, 80, 85, 86, 93, 94, 95, 96, 97, 98, 99, 100, 101], "Longest conjugated bond length": 58, "repeat_unit_match_count": 3, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}, {"Tag": 19442, "Chromophore": "CC(C)[Si](C#Cc%11c1cc(C(C)(C)C)ccc1c(C#Cc9c2ccc(C(C)(C)C)cc2c(C#Cc7c3cc(C(C)(C)C)ccc3c(C#Cc5c4ccc(C(C)(C)C)cc4c(C#C[Si](C(C)C)(C(C)C)C(C)C)c6ccc(C(C)(C)C)cc56)c8cc(C(C)(C)C)ccc78)c%10ccc(C(C)(C)C)cc9%10)c%12cc(C(C)(C)C)ccc%11%12)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 582.0, "Emission max (nm)": 609.0, "Lifetime (ns)": NaN, "Quantum yield": 0.02, "log(e/mol-1 dm3 cm-1)": 4.975431809, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 1139.68356, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.1302737505260163, "Emission max (eV)": 2.035828116266242, "Stokes shift (eV)": 0.09444563425977437, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96], "Longest conjugated bond length": 77, "repeat_unit_match_count": 4, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}, {"Tag": 19443, "Chromophore": "CC(C)[Si](C#Cc%15c2ccc(C(C)(C)C)cc2c(C#Cc%13c3cc(C(C)(C)C)ccc3c(C#Cc%11c4ccc(C(C)(C)C)cc4c(C#Cc9c5cc(C(C)(C)C)ccc5c(C#Cc7c6ccc(C(C)(C)C)cc6c(C#C[Si](C(C)C)(C(C)C)C(C)C)c8ccc(C(C)(C)C)cc78)c%10cc(C(C)(C)C)ccc9%10)c%12ccc(C(C)(C)C)cc%11%12)c%14cc(C(C)(C)C)ccc%13%14)c%16ccc(C(C)(C)C)cc%15%16)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 585.0, "Emission max (nm)": 623.0, "Lifetime (ns)": NaN, "Quantum yield": 0.018, "log(e/mol-1 dm3 cm-1)": 5.130333768, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 1339.92308, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.1193492697540877, "Emission max (eV)": 1.9900791698332927, "Stokes shift (eV)": 0.12927009992079497, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115], "Longest conjugated bond length": 96, "repeat_unit_match_count": 5, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}, {"Tag": 19444, "Chromophore": "CC(C)[Si](C#Cc%17c1cc(C(C)(C)C)ccc1c(C#Cc%15c2ccc(C(C)(C)C)cc2c(C#Cc%13c3cc(C(C)(C)C)ccc3c(C#Cc%11c4ccc(C(C)(C)C)cc4c(C#Cc9c5cc(C(C)(C)C)ccc5c(C#Cc7c6ccc(C(C)(C)C)cc6c(C#C[Si](C(C)C)(C(C)C)C(C)C)c8ccc(C(C)(C)C)cc78)c%10cc(C(C)(C)C)ccc9%10)c%12ccc(C(C)(C)C)cc%11%12)c%14cc(C(C)(C)C)ccc%13%14)c%16ccc(C(C)(C)C)cc%15%16)c%18cc(C(C)(C)C)ccc%17%18)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 589.0, "Emission max (nm)": 629.0, "Lifetime (ns)": NaN, "Quantum yield": NaN, "log(e/mol-1 dm3 cm-1)": 5.178976947, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 1540.1626, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.104956405443364, "Emission max (eV)": 1.9710959027124664, "Stokes shift (eV)": 0.13386050273089745, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134], "Longest conjugated bond length": 115, "repeat_unit_match_count": 6, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}], "data-23d98a78cf3a4d67e6e8e922ac3376db": [{"Type": "Emission max (eV)", "Shape": "circle"}, {"Type": "Absorption max (eV)", "Shape": "triangle"}]}}, {"mode": "vega-lite"});
</script>


For this series with 1-6 repeat units, the absorption (and emission) wavelength maxima indeed show a monotonic increase with increasing number of repeat units, which correlates to the length of the conjugated bond chain. So having a consistent molecular structure demonstrates that increasing the conjugated bond chain length increases the absorption maximum. Nagaoka et al. state this as "The bathochromic shifts [to higher wavelengths] in the UV−vis spectra suggested that the π-conjugation was extended with elongation of the linear chain." They thus experimentally investigated the question of whether the conjugated bond chain was extended by adding more repeat units, and concluded that it was. This is also expected from their
- X-ray structure of the `n = 3` oligomer: "The three anthracene plans are approximately coplanar along the linear molecular axis"
- Density Functional Theory (DFT) structure of all six oligomers: "the anthracene groups are coplanar...In these coplanar conformations, the molecular orbitals spread over the acetylene [adjacent carbons joined by a triple bond] and anthracene moieties at the HOMO and LUMO levels".

If the anthracene groups were not coplanar, the p orbitals on the acetylene and anthracene units would not be pointing in the same direction, potentially disrupting the conjugated bond chain. The RDKit method `bond.GetIsConjugated()` examines the molecular graph (with atoms as nodes and bonds as edges), which does not consider things like bond angles, so the experimental data more conclusively demonstrate that the conjugated bond chain length is increased with more repeat units.

Nagaoka et al. note "These results suggest that diacetylene linkers effectively mediate electronic communication between aromatic units".

As we discussed earlier, for comparisons such as the Stokes shift, it's better to plot in terms of energy. Let's do that and add the Stokes shift on the same x-axis, but below the main plot.


```python
# Order of curves from top to bottom
shapes = [absorption_shape, emission_shape, stokes_shape]

# Larger data point size for better visibility
point_size = 100

# First scatter plot for "Absorption max (eV)" and "Emission max (eV)" (upper part of the y-axis)
scatter_absorption = (
    alt.Chart(df_matches)
    .mark_point(size=point_size)
    .encode(
        x="repeat_unit_match_count",
        y=alt.Y("Absorption max (eV)", scale=alt.Scale(domain=[y_min_eV, y_max_eV])),
        color=alt.value("black"),
        shape=alt.value(absorption_shape),
        tooltip=["repeat_unit_match_count", "Absorption max (eV)"],
    )
)

scatter_emission = (
    alt.Chart(df_matches)
    .mark_point(size=point_size)
    .encode(
        x="repeat_unit_match_count",
        y=alt.Y(
            "Emission max (eV)", scale=alt.Scale(domain=[y_min_eV, y_max_eV])
        ),  # Same y-axis range as absorption
        color=alt.value("black"),
        shape=alt.value(emission_shape),
        tooltip=["repeat_unit_match_count", "Emission max (eV)"],
    )
)

# Second scatter plot for "Stokes shift (eV)" (lower part of the y-axis)
scatter_stokes = (
    alt.Chart(df_matches)
    .mark_point(size=point_size)
    .encode(
        x="repeat_unit_match_count",
        y=alt.Y(
            "Stokes shift (eV)",
            # Stokes shift on a different y-axis range
            scale=alt.Scale(
                domain=[0, stokes_y_max],
            ),
        ),
        color=alt.value("black"),
        shape=alt.value(stokes_shape),
        tooltip=["repeat_unit_match_count", "Stokes shift (eV)"],
    )
)

# Update legend data with shape information
legend_data = pl.DataFrame(
    {
        "Type": ["Absorption max (eV)", "Emission max (eV)", "Stokes shift (eV)"],
        "Shape": shapes,
    }
)

# Manual legend with explicit sort order
legend = (
    alt.Chart(legend_data)
    .mark_point(size=100)
    .encode(
        y=alt.Y(
            "Type:N",
            # Explicit sort order
            sort=[
                "Absorption max (eV)",
                "Emission max (eV)",
                "Stokes shift (eV)",
            ],
            axis=alt.Axis(orient="right", title=None),
        ),
        color=alt.value("black"),
        shape=alt.Shape(
            # Map the Shape column to the shapes in the legend
            "Shape:N",
            scale=alt.Scale(
                domain=shapes,
                range=shapes,
            ),
            # Disable shape legend
            legend=None,
        ),
    )
)

# Layer the absorption and emission scatter plots
upper_chart = alt.layer(scatter_absorption, scatter_emission).encode(
    y=alt.Y(
        "Absorption max (eV)",
        scale=alt.Scale(domain=[y_min_eV, y_max_eV]),
        title="Absorption or emission maximum (eV)",
    ),
    x=alt.X(
        "repeat_unit_match_count",
        scale=alt.Scale(domain=[0.5, max_repeat_count]),
        title=None,
        axis=alt.Axis(labels=False, tickSize=0),
    ),
)

# Create the Stokes shift scatter plot
lower_chart = scatter_stokes.encode(
    y=alt.Y(
        "Stokes shift (eV)",
        scale=alt.Scale(domain=[0, stokes_y_max]),
        title="Stokes shift (eV)",
    ),
    x=alt.X(
        "repeat_unit_match_count",
        scale=alt.Scale(domain=[0.5, max_repeat_count]),
        title="Number of repeat units",
    ),
)

# Combine the two charts vertically with some padding to simulate a break
combined_chart = alt.vconcat(
    upper_chart.properties(height=300), lower_chart.properties(height=100), spacing=10
).resolve_scale(x="shared")

# Display the combined chart with the legend
layered_chart = alt.hconcat(combined_chart, legend)

layered_chart.show()
```



<style>
  #altair-viz-702d7c6eba91438cbeebe0cbd66bc299.vega-embed {
    width: 100%;
    display: flex;
  }

  #altair-viz-702d7c6eba91438cbeebe0cbd66bc299.vega-embed details,
  #altair-viz-702d7c6eba91438cbeebe0cbd66bc299.vega-embed details summary {
    position: relative;
  }
</style>
<div id="altair-viz-702d7c6eba91438cbeebe0cbd66bc299"></div>
<script type="text/javascript">
  var VEGA_DEBUG = (typeof VEGA_DEBUG == "undefined") ? {} : VEGA_DEBUG;
  (function(spec, embedOpt){
    let outputDiv = document.currentScript.previousElementSibling;
    if (outputDiv.id !== "altair-viz-702d7c6eba91438cbeebe0cbd66bc299") {
      outputDiv = document.getElementById("altair-viz-702d7c6eba91438cbeebe0cbd66bc299");
    }
    const paths = {
      "vega": "https://cdn.jsdelivr.net/npm/vega@5?noext",
      "vega-lib": "https://cdn.jsdelivr.net/npm/vega-lib?noext",
      "vega-lite": "https://cdn.jsdelivr.net/npm/vega-lite@5.20.1?noext",
      "vega-embed": "https://cdn.jsdelivr.net/npm/vega-embed@6?noext",
    };

    function maybeLoadScript(lib, version) {
      var key = `${lib.replace("-", "")}_version`;
      return (VEGA_DEBUG[key] == version) ?
        Promise.resolve(paths[lib]) :
        new Promise(function(resolve, reject) {
          var s = document.createElement('script');
          document.getElementsByTagName("head")[0].appendChild(s);
          s.async = true;
          s.onload = () => {
            VEGA_DEBUG[key] = version;
            return resolve(paths[lib]);
          };
          s.onerror = () => reject(`Error loading script: ${paths[lib]}`);
          s.src = paths[lib];
        });
    }

    function showError(err) {
      outputDiv.innerHTML = `<div class="error" style="color:red;">${err}</div>`;
      throw err;
    }

    function displayChart(vegaEmbed) {
      vegaEmbed(outputDiv, spec, embedOpt)
        .catch(err => showError(`Javascript Error: ${err.message}<br>This usually means there's a typo in your chart specification. See the javascript console for the full traceback.`));
    }

    if(typeof define === "function" && define.amd) {
      requirejs.config({paths});
      require(["vega-embed"], displayChart, err => showError(`Error loading script: ${err.message}`));
    } else {
      maybeLoadScript("vega", "5")
        .then(() => maybeLoadScript("vega-lite", "5.20.1"))
        .then(() => maybeLoadScript("vega-embed", "6"))
        .catch(showError)
        .then(() => displayChart(vegaEmbed));
    }
  })({"config": {"view": {"continuousWidth": 300, "continuousHeight": 300}}, "hconcat": [{"vconcat": [{"layer": [{"mark": {"type": "point", "size": 100}, "encoding": {"color": {"value": "black"}, "shape": {"value": "triangle"}, "tooltip": [{"field": "repeat_unit_match_count", "type": "quantitative"}, {"field": "Absorption max (eV)", "type": "quantitative"}], "x": {"field": "repeat_unit_match_count", "type": "quantitative"}, "y": {"field": "Absorption max (eV)", "scale": {"domain": [1.872541107576843, 2.958659747605565]}, "type": "quantitative"}}}, {"mark": {"type": "point", "size": 100}, "encoding": {"color": {"value": "black"}, "shape": {"value": "circle"}, "tooltip": [{"field": "repeat_unit_match_count", "type": "quantitative"}, {"field": "Emission max (eV)", "type": "quantitative"}], "x": {"field": "repeat_unit_match_count", "type": "quantitative"}, "y": {"field": "Emission max (eV)", "scale": {"domain": [1.872541107576843, 2.958659747605565]}, "type": "quantitative"}}}], "encoding": {"x": {"axis": {"labels": false, "tickSize": 0}, "field": "repeat_unit_match_count", "scale": {"domain": [0.5, 6.5]}, "title": null, "type": "quantitative"}, "y": {"field": "Absorption max (eV)", "scale": {"domain": [1.872541107576843, 2.958659747605565]}, "title": "Absorption or emission maximum (eV)", "type": "quantitative"}}, "height": 300}, {"mark": {"type": "point", "size": 100}, "encoding": {"color": {"value": "black"}, "shape": {"value": "square"}, "tooltip": [{"field": "repeat_unit_match_count", "type": "quantitative"}, {"field": "Stokes shift (eV)", "type": "quantitative"}], "x": {"field": "repeat_unit_match_count", "scale": {"domain": [0.5, 6.5]}, "title": "Number of repeat units", "type": "quantitative"}, "y": {"field": "Stokes shift (eV)", "scale": {"domain": [0, 0.1472465530039872]}, "title": "Stokes shift (eV)", "type": "quantitative"}}, "height": 100}], "data": {"name": "data-fcb24888a4ab58ce1808fb0997ef05d4"}, "resolve": {"scale": {"x": "shared"}}, "spacing": 10}, {"data": {"name": "data-b6fc9c0ba0914d2821e57e91562567a2"}, "mark": {"type": "point", "size": 100}, "encoding": {"color": {"value": "black"}, "shape": {"field": "Shape", "legend": null, "scale": {"domain": ["triangle", "circle", "square"], "range": ["triangle", "circle", "square"]}, "type": "nominal"}, "y": {"axis": {"orient": "right", "title": null}, "field": "Type", "sort": ["Absorption max (eV)", "Emission max (eV)", "Stokes shift (eV)"], "type": "nominal"}}}], "$schema": "https://vega.github.io/schema/vega-lite/v5.20.1.json", "datasets": {"data-fcb24888a4ab58ce1808fb0997ef05d4": [{"Tag": 19439, "Chromophore": "CC(C)[Si](C#Cc1c2ccc(C(C)(C)C)cc2c(C#C[Si](C(C)C)(C(C)C)C(C)C)c2ccc(C(C)(C)C)cc12)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 440.0, "Emission max (nm)": 461.0, "Lifetime (ns)": NaN, "Quantum yield": 0.92, "log(e/mol-1 dm3 cm-1)": 4.413299764, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 651.18004, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.817771188195776, "Emission max (eV)": 2.6894128477356647, "Stokes shift (eV)": 0.1283583404601112, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 29, 30, 31, 32, 37, 38, 45, 46, 47], "Longest conjugated bond length": 20, "repeat_unit_match_count": 1, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}, {"Tag": 19440, "Chromophore": "CC(C)[Si](C#Cc5c1cc(C(C)(C)C)ccc1c(C#Cc3c2ccc(C(C)(C)C)cc2c(C#C[Si](C(C)C)(C(C)C)C(C)C)c4ccc(C(C)(C)C)cc34)c6cc(C(C)(C)C)ccc56)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 523.0, "Emission max (nm)": 541.0, "Lifetime (ns)": NaN, "Quantum yield": 0.02, "log(e/mol-1 dm3 cm-1)": 4.694605199, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 963.6346, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.3705914393998877, "Emission max (eV)": 2.29171778707235, "Stokes shift (eV)": 0.07887365232753751, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 28, 29, 30, 31, 32, 43, 44, 45, 46, 51, 52, 53, 54, 55, 56, 61, 62, 69, 70, 71, 72, 73, 74], "Longest conjugated bond length": 39, "repeat_unit_match_count": 2, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}, {"Tag": 19441, "Chromophore": "CC(C)[Si](C#Cc8c1cc(C(C)(C)C)ccc1c(C#Cc6c2ccc(C(C)(C)C)cc2c(C#Cc4c3cc(C(C)(C)C)ccc3c(C#C[Si](C(C)C)(C(C)C)C(C)C)c5cc(C(C)(C)C)ccc45)c7ccc(C(C)(C)C)cc67)c9cc(C(C)(C)C)ccc89)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 570.0, "Emission max (nm)": 589.0, "Lifetime (ns)": NaN, "Quantum yield": 0.022, "log(e/mol-1 dm3 cm-1)": 4.781036939, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 1276.08916, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.175121618958143, "Emission max (eV)": 2.104956405443364, "Stokes shift (eV)": 0.070165213514779, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 42, 43, 44, 45, 46, 57, 58, 59, 60, 65, 66, 67, 68, 69, 70, 75, 76, 77, 78, 79, 80, 85, 86, 93, 94, 95, 96, 97, 98, 99, 100, 101], "Longest conjugated bond length": 58, "repeat_unit_match_count": 3, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}, {"Tag": 19442, "Chromophore": "CC(C)[Si](C#Cc%11c1cc(C(C)(C)C)ccc1c(C#Cc9c2ccc(C(C)(C)C)cc2c(C#Cc7c3cc(C(C)(C)C)ccc3c(C#Cc5c4ccc(C(C)(C)C)cc4c(C#C[Si](C(C)C)(C(C)C)C(C)C)c6ccc(C(C)(C)C)cc56)c8cc(C(C)(C)C)ccc78)c%10ccc(C(C)(C)C)cc9%10)c%12cc(C(C)(C)C)ccc%11%12)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 582.0, "Emission max (nm)": 609.0, "Lifetime (ns)": NaN, "Quantum yield": 0.02, "log(e/mol-1 dm3 cm-1)": 4.975431809, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 1139.68356, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.1302737505260163, "Emission max (eV)": 2.035828116266242, "Stokes shift (eV)": 0.09444563425977437, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96], "Longest conjugated bond length": 77, "repeat_unit_match_count": 4, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}, {"Tag": 19443, "Chromophore": "CC(C)[Si](C#Cc%15c2ccc(C(C)(C)C)cc2c(C#Cc%13c3cc(C(C)(C)C)ccc3c(C#Cc%11c4ccc(C(C)(C)C)cc4c(C#Cc9c5cc(C(C)(C)C)ccc5c(C#Cc7c6ccc(C(C)(C)C)cc6c(C#C[Si](C(C)C)(C(C)C)C(C)C)c8ccc(C(C)(C)C)cc78)c%10cc(C(C)(C)C)ccc9%10)c%12ccc(C(C)(C)C)cc%11%12)c%14cc(C(C)(C)C)ccc%13%14)c%16ccc(C(C)(C)C)cc%15%16)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 585.0, "Emission max (nm)": 623.0, "Lifetime (ns)": NaN, "Quantum yield": 0.018, "log(e/mol-1 dm3 cm-1)": 5.130333768, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 1339.92308, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.1193492697540877, "Emission max (eV)": 1.9900791698332927, "Stokes shift (eV)": 0.12927009992079497, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115], "Longest conjugated bond length": 96, "repeat_unit_match_count": 5, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}, {"Tag": 19444, "Chromophore": "CC(C)[Si](C#Cc%17c1cc(C(C)(C)C)ccc1c(C#Cc%15c2ccc(C(C)(C)C)cc2c(C#Cc%13c3cc(C(C)(C)C)ccc3c(C#Cc%11c4ccc(C(C)(C)C)cc4c(C#Cc9c5cc(C(C)(C)C)ccc5c(C#Cc7c6ccc(C(C)(C)C)cc6c(C#C[Si](C(C)C)(C(C)C)C(C)C)c8ccc(C(C)(C)C)cc78)c%10cc(C(C)(C)C)ccc9%10)c%12ccc(C(C)(C)C)cc%11%12)c%14cc(C(C)(C)C)ccc%13%14)c%16ccc(C(C)(C)C)cc%15%16)c%18cc(C(C)(C)C)ccc%17%18)(C(C)C)C(C)C", "Solvent": "Cc1ccccc1", "Absorption max (nm)": 589.0, "Emission max (nm)": 629.0, "Lifetime (ns)": NaN, "Quantum yield": NaN, "log(e/mol-1 dm3 cm-1)": 5.178976947, "abs FWHM (cm-1)": NaN, "emi FWHM (cm-1)": NaN, "abs FWHM (nm)": NaN, "emi FWHM (nm)": NaN, "Molecular weight (g mol-1)": 1540.1626, "Reference": "DOI: 10.1021/acs.joc.8b00311", "Absorption max (eV)": 2.104956405443364, "Emission max (eV)": 1.9710959027124664, "Stokes shift (eV)": 0.13386050273089745, "longest_bond_indices": [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134], "Longest conjugated bond length": 115, "repeat_unit_match_count": 6, "terminal_match_count": 2, "too_many_rings_match_count": 0, "too_many_rings_N_match_count": 0}], "data-b6fc9c0ba0914d2821e57e91562567a2": [{"Type": "Absorption max (eV)", "Shape": "triangle"}, {"Type": "Emission max (eV)", "Shape": "circle"}, {"Type": "Stokes shift (eV)", "Shape": "square"}]}}, {"mode": "vega-lite"});
</script>


The y-axis makes it easier to compare energy values both between and within molecules. This plot agrees with Nagaoka et al's statement that "HOMO—LUMO gaps [absorption maxima, in eV] decrease with an increasing chain length".

The Stokes shift is relatively constant, ranging between 0.7 and 0.14 eV, despite the absorption maximum decreasing by ~0.7 eV from `n = 1` to `n = 6`, because the emission maximum follows the same trend. The Stokes shift does show a trend, decreasing from `n = 1` until `n = 3`, then increasing up to `n = 6`, which may be due to conflicting factors.

## Conclusions

By finding a series of molecules with a consistent structure, and from 1-6 repeat units, we found that absorption maximum increases (from the ultraviolet to the visible) as the length of the conjugated bond chain increases, as predicted by molecular orbital theory. This effect helps explain why some organic molecules are colored.

By changing the conjugated bond chain length, chemists can tune the optical properties of absorption and emission wavelength to make useful devices. It's important to verify experimentally that the absorption wavelength maximum does in fact increase, indicating that the electrons can spread over more bonds, when the chain length (as determined by cheminformatics) is increased.