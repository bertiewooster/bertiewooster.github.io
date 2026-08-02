import marimo

__generated_with = "0.23.9"
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

    return


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
    This is a simple one, but helpful if you want to let the user tab between different outputs, in this case two entity-relationship diagrams created by graphviz Digraph:
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ![Entity-relationship diagrams in Marimo tabs](../images/ChEMBL-Marimo-post/ChEMBL-Marimo-tabs.png)
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Selecting Multiple Items From a Menu
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The lets you select multiple chemical compounds for lookup in ChEMBL.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ![Drop-down multiselect menu for ChEMBL compounds](../images/ChEMBL-Marimo-post/ChEMBL-Marimo-dropdown-multiselect.png)
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Action Button
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Clicking this button make the call to ChEBML API via its chembl_webresource_client to fetch data for the selected compounds.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ![Button Fetch from ChEMBL](../images/ChEMBL-Marimo-post/ChEMBL-Marimo-fetch-from-ChEMBL.png)
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
 
    """)
    return


if __name__ == "__main__":
    app.run()
