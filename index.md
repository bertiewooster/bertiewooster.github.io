Scientific software developer in the Washington, D.C. area.

# Portfolio of my projects

## [RDKit Utility to Check Whether Starting Materials for Synthesizing Your Target Molecules Are Commercially Available]({% post_url 2023-02-07-Are-the-Starting-Materials-for-Synthesizing-Your-Target-Molecules-Commercially-Available %})
*Uses Python, RDKit, PubChem's API, asyncio, and Semaphore*

[<img alt="Three reactions, each in a row. First column: Target molecule and whether it's accessible based on commercial availability of reactants. Subsequent columns: Each reactant and whether it's commercial available." style="width:600px; height:600px" src="/images/reaction-accessible.png">]({% post_url 2023-02-07-Are-the-Starting-Materials-for-Synthesizing-Your-Target-Molecules-Commercially-Available %})

Given target molecules and reactions to synthesize them, determine whether the starting materials are commercially available using PubChem's API, and thus whether the target is synthetically accessible.

## [RDKit Utility to Create a Mass Spectrometry Fragmentation Tree]({% post_url 2023-01-02-Mass-Spectrometry-Fragmentation-Tree %})
*Uses Python and RDKit*

[<img alt="Annotated mass spectrometry fragmentation tree using the function mass_spec_frag_tree in this blog post" style="width:518px; height:393px" src="/images/mass_spec_frag_tree_CAIr_annotated.png">]({% post_url 2023-01-02-Mass-Spectrometry-Fragmentation-Tree %})

Given a mass spec fragmentation hierarchy, with species as SMILES strings, display the fragmentation tree in a grid, labeling each species with its name and either mass or mass to charge ratio `m/z`.

## [RDKit Utility to Find the Maximum Common Substructure, and Groups Off It, Between a Set of Molecules]({% post_url 2022-12-25-RDKit-Find-Groups-Off-Common-Core %})
*Uses Python and RDKit*

[<img alt="Annotated grid of maximum common substructure and core; molecules and groups off maximum common substructure" style="width:600px; height:551px" src="/images/pyridinols-MCS-groups-annotated.png">]({% post_url 2022-12-25-RDKit-Find-Groups-Off-Common-Core %})

Given a collection of molecules as SMILES strings, find the maximum common substructure (MCS) match between them, and the groups off that common core for each molecule, displaying the results using a grid.

## [Chemistry machine learning for drug discovery with DeepChem]({% post_url 2022-12-13-Chemistry-machine-learning-for-drug-discovery-with-DeepChem %})
*Uses Python, DeepChem, seaborn, Matplotlib, and pandas*

[<img alt="Predicted against measured lipophilicty for test and train data" style="width:384px; height:262px" src="/images/lipophilicity-machine-learning-graph.png">]({% post_url 2022-12-13-Chemistry-machine-learning-for-drug-discovery-with-DeepChem %})

Use the DeepChem deep learning package to predict compounds' lipophilicity--how well they are absorbed into the lipids of biological membranes, which is important for oral delivery of drugs.

## [RDKit Utility to Visualize Retrosynthetic Analysis Hierarchically]({% post_url 2022-11-11-RDKit-Recap-decomposition-tree %})
*Uses Python and RDKit*

[<img alt="Annotated Recap retrosynthetic hierarchy tree" style="width:513px; height:418px" src="/images/tree_marked_up.png">]({% post_url 2022-10-09-RDKit-find-and-highlight-the-maximum-common-substructure-between-molecules %})

Given a target molecule, use the [Recap algorithm](https://www.semanticscholar.org/paper/RECAP-%E2%80%94-Retrosynthetic-Combinatorial-Analysis-A-New-Lewell-Judd/fbfb10d1f63aa803f6d47df6587aa0e41109f5ee){:target='_blank'} to decompose it into a set of fragments that could be combined to make the parent molecule using common reactions. Display the fragmentation hierarchically.

## [RDKit Utility to Find and Highlight the Maximum Common Substructure Amongst Molecules]({% post_url 2022-10-09-RDKit-find-and-highlight-the-maximum-common-substructure-between-molecules %})
*Uses Python and RDKit*

[<img alt="Maximum substructure match, and the two molecules which are labeled by their functional groups" style="width:600px; height:200px" src="/images/RDKit-MCS-example.png">]({% post_url 2022-10-09-RDKit-find-and-highlight-the-maximum-common-substructure-between-molecules %})

Given a collection of molecules as [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system){:target='_blank'} strings, find the maximum common substructure (MCS) match between them as a [SMARTS](https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification){:target='_blank'} string, display the match pattern as a molecule, and highlight the match pattern in each molecule using a grid.

## [Materials and Cheminformatics Sampler](https://sampler-flask.herokuapp.com/)
*Uses Python, NumPy, SymPy, ChemPy, Flask, JavaScript, and Bootstrap*

[<img style="width:400px; height:229px;" src="images/sampler-constraints.png">](https://sampler-flask.herokuapp.com/)

Find a given number of points which satisfy constraints given in a constraints file for an n-dimensional space defined on the unit hypercube, then write them to an output file.

Optionally, identify the components (dimensions) in the constraints file using chemical formulas, and Sampler will use ChemPy to calculate their molar masses, then output the component weight fraction.

[<img style="float: right; width:130px; height:85px;" src="images/ptable-terms-highlighted.png">](https://ptablenav.herokuapp.com/)
## [Periodic Table Navigator](https://ptablenav.herokuapp.com/)
*Uses Ruby, Sinatra, PostgreSQL, and JavaScript*

Understand how the elements are related to each other. Emphasizes electronic configuration of the elements.

# My open-source contributions

## RDKit cheminformatics package
- [Improved documentation](https://github.com/rdkit/rdkit/pulls?q=is%3Amerged+is%3Apr+author%3Abertiewooster+) by illustrating drawing capability in tutorial and adding SMILES (chemical notation) for R groups

## SymPy computer algebra system in pure Python
- Technical writer for funded 2022 Season of Docs project: Creating documentation for how to [solve equations](https://docs.sympy.org/dev/guides/solving/index.html)
- Core developer wrote "I think you are doing excellent work on the SymPy documentation. Thank you!"
- [Led selection of new Sphinx theme for SymPy documentation](https://github.com/sympy/sympy/issues/22716); the [new theme was implemented](https://docs.sympy.org/dev/)
- [Contributed code for documentation](https://github.com/sympy/sympy/pulls?q=is:pr+author:bertiewooster+is:merged) to explain usage of a core class for users and developers, and improve accessibility
- Lead developer wrote “You've been doing great work with the Sphinx theme and other documentation work”

## ChemPy package for chemistry in Python
- Initiated and provided scientific and coding direction to issue to [improve interpretation of chemical formulas](https://github.com/bjodah/chempy/issues/202)
- Spurred a developer to improve code
- Package author wrote “Great work guys!”

## Sphinx documentation generator
- Initiated issue to [improve accessibility and internationalization of documentation generated by Sphinx](https://github.com/sphinx-doc/sphinx/issues?q=author%3Abertiewooster+); was addressed within a day by Sphinx’s main developer
