import marimo

__generated_with = "0.18.4"
app = marimo.App()


@app.cell
def _():
    import matplotlib.pyplot as plt
    return (plt,)


@app.cell
def _():
    a = 7
    return (a,)


@app.cell
def _(a):
    b = a + 3
    b
    return (b,)


@app.cell
def _(a, b, plt):
    # Create the plot
    fig, ax = plt.subplots()
    ax.scatter(a, b)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_title("Plot of b against a")

    fig
    return


if __name__ == "__main__":
    app.run()
