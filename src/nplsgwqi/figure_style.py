from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.ticker import MaxNLocator

NATURE_RCPARAMS = {
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "axes.linewidth": 0.8,
    "lines.linewidth": 2.0,
    "lines.markersize": 6,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    "grid.linestyle": "--",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
}


def apply_publication_style() -> None:
    mpl.rcParams.update(NATURE_RCPARAMS)


def style_axes(
    ax: Axes,
    *,
    xbins: int | None = None,
    ybins: int | None = None,
    grid: bool = False,
    grid_axis: str = "both",
) -> Axes:
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.tick_params(which="both", direction="in", top=True, right=True, length=6, width=0.8)
    ax.tick_params(which="minor", length=3, width=0.6)
    if xbins is not None:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=xbins))
    if ybins is not None:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=ybins))
    ax.grid(False)
    if grid:
        ax.grid(True, axis=grid_axis, alpha=0.3, linewidth=0.5, linestyle="--")
    return ax


def style_colorbar(cbar, label: str | None = None) -> None:
    if label is not None:
        cbar.set_label(label, fontsize=16)
    cbar.ax.tick_params(which="both", direction="in", labelsize=14)


def style_legend(legend: Legend | None) -> Legend | None:
    if legend is None:
        return None
    frame = legend.get_frame()
    frame.set_alpha(0.8)
    frame.set_edgecolor("0.8")
    frame.set_linewidth(0.8)
    return legend


def save_figure(path: Path, *, fig: Figure | None = None, dpi: int = 600, pad: float = 0.35) -> None:
    fig = fig or plt.gcf()
    path.parent.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*not compatible with tight_layout.*", category=UserWarning)
        fig.tight_layout(pad=pad)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
