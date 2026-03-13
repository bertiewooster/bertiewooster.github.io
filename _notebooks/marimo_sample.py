import marimo

__generated_with = "0.19.11"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return


@app.cell
def _():
    a = 7
    return (a,)


@app.cell
def _(a):
    a + 2
    return


if __name__ == "__main__":
    app.run()
