import marimo

__generated_with = "0.23.6"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Interactive Marimo Notebook for Drug-Like ChEMBL Compounds Within Target Profiles
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    I published [Prioritizing Drug-Like ChEMBL Compounds Within Target Profiles]({% post_url 2026-01-05-ChEBML-database %}) as my first Marimo notebook. I largely used it as a replacement for a Jupyter Notebook, which my previous blog posts had been.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    After I published the blog post, the Marimo team contacted me about presenting at their Community Call. They liked how I approached the virtual screening workflow by looking into target profiles.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The Marimo team proposed adding my notebook to **mo**lab, Marimo's online platform for running notebooks in a web browser. At first I didn't think this was going to work because [pyodide doesn't yet support RDKit](https://github.com/pyodide/pyodide-recipes/issues/512) and I was trying to run Molab in WASM/Pyodide execution mode. But Marimo told me about the other execution mode, a cloud backend that runs full Python with any package, and that worked.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The Marimo team made a more interactive version of the notebook, then I made some changes. The result: [🧪 Prioritizing Drug-Like ChEMBL Compounds Within Target Profiles](https://molab.marimo.io/notebooks/nb_pFcNtsNCS1Qy2Ua4riPzqK)
    """)
    return


@app.cell
def _():
    from graphviz import Digraph

    return (Digraph,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Interactive Elements
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Tabs
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    This is a simple one, but helpful if you want to let the user tab between different outputs. Here we use graphviz Digraph to create two diagrams:
    """)
    return


@app.cell
def _(Digraph):
    _kw = dict(fontsize="8", arrowhead="crow", arrowtail="crow",
               dir="both", color="black", arrowsize="0.4")

    dot_chembl = Digraph(format="svg")
    dot_chembl.attr(rankdir="LR", splines="ortho")
    dot_chembl.attr("node", shape="record", fontsize="10", style="filled")
    dot_chembl.node("Compound", "{Compound|compound_id (PK)}", fillcolor="#A3C1DA")
    dot_chembl.node("Activity", "{Activity|activity_id (PK)}", fillcolor="#A3C1DA")
    dot_chembl.node("Target",   "{Target|target_id (PK)}",    fillcolor="#A3C1DA")
    dot_chembl.edge("Compound", "Activity", **_kw)
    dot_chembl.edge("Activity", "Target",   **_kw)

    dot_simple = Digraph(format="svg")
    dot_simple.attr(rankdir="LR", splines="ortho")
    dot_simple.attr("node", shape="record", fontsize="10", style="filled")
    dot_simple.node("Compound", "{Compound|compound_id (PK)}", fillcolor="#A3C1DA")
    dot_simple.node("Target",   "{Target|target_id (PK)}",    fillcolor="#A3C1DA")
    dot_simple.edge("Compound", "Target", **_kw)
    return dot_chembl, dot_simple


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    then set up Marimo tabs to let the user click a tab title to go to that diagram:
    """)
    return


@app.cell
def _(dot_chembl, dot_simple, mo):
    mo.ui.tabs({
        "ChEMBL: Compound ↔ Activity ↔ Target": mo.Html(dot_chembl.pipe(format="svg").decode()),
        "Our schema: Compound ↔ Target (direct)": mo.Html(dot_simple.pipe(format="svg").decode()),
    })
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Selecting Multiple Items From a Menu
    """)
    return


@app.cell
def _(mo):
    _curated_ids = [
        "CHEMBL25",       # Aspirin
        "CHEMBL112",      # Acetaminophen
        "CHEMBL941",      # Imatinib (Gleevec)
        "CHEMBL1737",     # Sildenafil (Viagra)
        "CHEMBL185",      # Fluorouracil
        "CHEMBL41",       # Fluoxetine (Prozac)
        "CHEMBL113",      # Caffeine
        "CHEMBL5416410",  # Dasatinib
        "CHEMBL1079742",  # Erlotinib (Tarceva)
        "CHEMBL939",      # Gefitinib (Iressa)
        "CHEMBL796",      # Methylphenidate
        "CHEMBL809",      # Setraline
    ]
    compound_ids_ui = mo.ui.multiselect(
        options=_curated_ids,
        value=["CHEMBL25", "CHEMBL796", "CHEMBL809"],
        label="Select compounds to fetch",
    )
    compound_ids_ui
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
