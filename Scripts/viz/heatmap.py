"""
Reusable heatmap helpers for the compartment-connectivity figures.
"""
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list

sns.set_style("ticks")

#------------------------------------------------------------------------------
def ward_order(matrix, method = "ward", metric = "euclidean"):
    """Return row order from Ward clustering of [matrix | matrix.T]."""
    fused        = pd.concat((matrix, matrix.T), axis = 1)
    link         = linkage(fused.fillna(0), method = method, metric = metric)
    return leaves_list(link)


#------------------------------------------------------------------------------
def plot_connectivity_heatmap(df,
                              norm,
                              ax,
                              cluster = True,
                              order   = None,
                              label   = True,
                              cmap    = "magma",
                              shrink  = 0.2):
    """
    Plot a square pre/post connectivity heatmap.

    Parameters
    ----------
    df : pandas.DataFrame
        Long-form with columns ``Pre``, ``Post``, ``weight`` (or
        ``Source_Compartment`` / ``Target_Compartment`` / ``Weight``).
    norm : matplotlib normalization (e.g. PowerNorm).
    ax : matplotlib axis.
    cluster : bool
        If True, Ward-cluster rows/columns. Ignored when ``order`` is given.
    order : list, optional
        Explicit row/column order (overrides clustering).
    label : bool
        Whether to show tick labels.
    """
    # Standardize column names
    df = df.rename(columns = {"Source_Compartment" : "Pre",
                              "Target_Compartment" : "Post",
                              "Weight"             : "weight"})

    items  = sorted(set(df["Pre"]) | set(df["Post"]))
    matrix = (df.pivot(index = "Pre", columns = "Post", values = "weight")
                .reindex(index = items, columns = items)
                .fillna(0))

    if order is not None:
        keep   = [o for o in order if o in matrix.index]
        matrix = matrix.reindex(index = keep, columns = keep).fillna(0)
    elif cluster:
        idx    = ward_order(matrix)
        matrix = matrix.iloc[idx, idx]

    sns.heatmap(matrix,
                cmap        = cmap,
                norm        = norm,
                ax          = ax,
                square      = True,
                cbar_kws    = {"label"  : "Connection Weight",
                               "shrink" : shrink,
                               "aspect" : 20},
                xticklabels = label,
                yticklabels = label)

    if label:
        ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, fontsize = 4)
        ax.set_yticklabels(ax.get_yticklabels(), rotation = 0,  fontsize = 4)

    return matrix
