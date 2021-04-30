"""
Regression plot (experimental vs predicted).
"""

import os
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
from sklearn.metrics import mean_squared_error

def regplot(
    df
) -> None:

    color = "orange"

    true = df["true"].to_numpy()
    predicted = df["predicted"].to_numpy()

    m = min(min(true), min(predicted))
    M = max(max(true), max(predicted))

    pr = stats.pearsonr(true, predicted)[0]
    rmse = np.sqrt(mean_squared_error(true, predicted))

    g = sns.jointplot(x=true, y=predicted, kind="reg", color=color)

    # Add std as errorbar
    g.ax_joint.errorbar(
        true, predicted, yerr=df["std"].to_numpy(), fmt="none", ecolor=color, alpha=0.25
    )

    x = np.linspace(m, M, 1000)
    g.ax_joint.plot(x, x, color="k", alpha=0.5, zorder=-10)
    g.ax_joint.fill_between(x, x + 1, x - 1, alpha=0.2, color="grey", zorder=-10)
    g.ax_joint.fill_between(x, x + 2, x - 2, alpha=0.1, color="grey", zorder=-10)

    g.ax_joint.set_xlabel("Experimental")
    g.ax_joint.set_ylabel("Predicted")

    g.ax_joint.text(
        0.5,
        0.05,
        f"RMSE = {rmse:.2f} | Pearson's $r$ = {pr:.2f}",
        horizontalalignment="center",
        verticalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.5},
        transform=g.ax_joint.transAxes,
    )

    plt.tight_layout()

    for ext in [".pdf", ".png"]:
        plt.savefig(f"images/regplot{ext}")

if __name__ == "__main__":

    # Create fake data
    x = np.linspace(-10, 10, 50)
    x += np.random.normal(0, 1, size=x.shape)
    y = {i: x + np.random.normal(0, 5, size=x.shape) for i in range(5)}

    # Store fake data with mean and standard deviation
    df = pd.DataFrame.from_dict(y)
    df["predicted"] = df.mean(axis="columns")
    df["std"] = df.drop(columns="predicted").std(axis="columns")
    df["true"] = x

    regplot(df)
