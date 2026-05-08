#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sensory (S-* and BolwigsAxon) lineage innervation density at level2.

Whole-neuron T-bar and PSD distributions across Ipsi/Contra-collapsed
compartments, row-normalized and ordered by fused Ward clustering.

Inputs (read-only):
    $PROJECTS_HOME/L1_Neuron-Typing/Analysis_Outputs/Synapses/
        Synapse-by-Lineage_level2.parquet

Outputs (Analysis_Outputs/Sensory_Innervation/):
    Sensory_Innervation_TBars_level2.parquet
    Sensory_Innervation_PSDs_level2.parquet
"""
import numpy as np
import pandas as pd
import os

from scipy.cluster.hierarchy import leaves_list, linkage

#------------------------------------------------------------------------------
FOLDER = "Analysis_Outputs/Sensory_Innervation"
HOME   = os.environ["PROJECTS_HOME"]

INCLUDE_PATTERNS = ("S-", "BolwigsAxon")

#------------------------------------------------------------------------------
def add_relative_compartment(df):
    """
    Collapse left/right hemispheres into Ipsi/Contra.

    """
    df = df.copy()

    bare = df["Compartment"].str.rsplit("_", n = 1).str[0]
    side = df["Compartment"].str.rsplit("_", n = 1).str[1]

    df["CompRelative"] = np.where(df["Hemisphere"] == side, "Ipsi", "Contra")
    df["Comp"]         = bare + "_" + df["CompRelative"]
    return df

#------------------------------------------------------------------------------
def build_matrix(df, lineages):
    """
    Pivot to (Lineage x Comp), average across hemispheres and row-normalize.

    """
    avg = (df.groupby(["Lineage", "Comp"])["Count"]
             .mean()
             .reset_index())

    matrix = avg.pivot(index = "Lineage", columns = "Comp", values = "Count")
    matrix = matrix.reindex(index = lineages).fillna(0)

    matrix = matrix.div(matrix.sum(axis = 1), axis = 0)
    return matrix

#------------------------------------------------------------------------------
def ward_reorder(*matrices):
    """
    Ward-cluster the row-fusion of the given matrices and apply the leaf
    order (reversed, matching the notebooks) to each one.

    """
    fused = pd.concat(matrices, axis = 1, join = "inner").fillna(0)
    order = leaves_list(linkage(fused, method = "ward", metric = "euclidean"))[::-1]

    return [m.loc[fused.index[order]] for m in matrices]

#------------------------------------------------------------------------------
if __name__ == "__main__":
    # Load level2 lineage table
    df = pd.read_parquet(f"{HOME}/L1_Neuron-Typing/Analysis_Outputs/Synapses/Synapse-by-Lineage_level2.parquet")

    # Keep only sensory and BolwigsAxon lineages
    df = df.loc[df["Lineage"].str.contains("|".join(INCLUDE_PATTERNS))]
    df = add_relative_compartment(df)

    lineages = np.unique(df["Lineage"])

    tbar_matrix = build_matrix(df.loc[df["Connector"] == "TBars"], lineages)
    psd_matrix  = build_matrix(df.loc[df["Connector"] == "PSDs"],  lineages)

    # Fused Ward clustering — same row order on both matrices
    tbar_matrix, psd_matrix = ward_reorder(tbar_matrix, psd_matrix)

    tbar_matrix.reset_index().to_parquet(f"{FOLDER}/Sensory_Innervation_TBars_level2.parquet", index = False, engine = "pyarrow")
    psd_matrix.reset_index().to_parquet(f"{FOLDER}/Sensory_Innervation_PSDs_level2.parquet",  index = False, engine = "pyarrow")

    print("हो गया दोस्तों!")
