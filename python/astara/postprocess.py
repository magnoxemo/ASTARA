"""Post-processing helpers for `Simulation.run()` output.

Thin, Jupyter-friendly wrappers around pandas/matplotlib -- nothing here
does anything a notebook couldn't do in a couple of lines itself, they just
save writing those lines every time. Requires the `astara[analysis]` extra
(`pip install astara[analysis]`); each function raises a clear `ImportError`
if it's missing rather than failing on an unrelated `NameError`.
"""

from __future__ import annotations

import importlib


def _require(module_name: str, package_hint: str):
    try:
        return importlib.import_module(module_name)
    except ImportError as exc:
        raise ImportError(
            f"{module_name} is required for astara.postprocess.{package_hint}; "
            "install it with `pip install astara[analysis]`."
        ) from exc


def plot_traces(df, columns, x: str = "t_s", title: str | None = None, ax=None, **plot_kwargs):
    """Line-plot `columns` of `df` against `x`. Returns the `Axes` used."""
    plt = _require("matplotlib.pyplot", "plot_traces")

    if ax is None:
        _, ax = plt.subplots()
    for column in columns:
        ax.plot(df[x], df[column], label=column, **plot_kwargs)
    ax.set_xlabel(x)
    ax.legend()
    if title:
        ax.set_title(title)
    return ax


def compare(baseline, other, columns, x: str = "t_s"):
    """Row-align two `Simulation.run()` DataFrames on `x` and return `other - baseline`
    for each column in `columns` (plus the shared `x` column) -- e.g. comparing a
    perturbed run against an untouched baseline run of the same scenario."""
    pd = _require("pandas", "compare")

    merged = pd.merge(
        baseline[[x] + list(columns)],
        other[[x] + list(columns)],
        on=x,
        suffixes=("_baseline", "_other"),
    )
    result = pd.DataFrame({x: merged[x]})
    for column in columns:
        result[column] = merged[f"{column}_other"] - merged[f"{column}_baseline"]
    return result


def to_csv(df, path, **kwargs) -> None:
    df.to_csv(path, index=False, **kwargs)


def to_parquet(df, path, **kwargs) -> None:
    df.to_parquet(path, **kwargs)
