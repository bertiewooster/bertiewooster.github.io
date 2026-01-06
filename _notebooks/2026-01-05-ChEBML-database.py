import marimo

__generated_with = "0.18.4"
app = marimo.App()


@app.cell
def _():
    a = 7
    return (a,)


@app.cell
def _(a):
    a + 3
    return


if __name__ == "__main__":
    app.run()
