#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lineage compartment-innervation matrices and focusness metrics — Shannon
entropy, Simpson (HHI) index, their effective-number-of-compartments forms,
and hierarchy / breadth ratios — computed per lineage from the synapse
tables produced by L1_Neuron-Typing.

For each level the script writes:

* Ipsi/Contra-collapsed (Lineage x Comp) matrices, row-normalized, ordered
  by Ward clustering (fused TBar+PSD linkage for the whole matrices,
  PSD-only linkage for the dendrite matrix).
* Per-lineage focusness metrics derived from those matrices.

The dendrite output uses PSDs only — TBars are dropped before any matrix
is built.

Inputs (read-only):
    $PROJECTS_HOME/L1_Neuron-Typing/Analysis_Outputs/Synapses/
        Synapse-by-Lineage_level1.parquet
        Synapse-by-Lineage_level2.parquet
        Synapse-by-Lineage_level1_Fragments.parquet
        Synapse-by-Lineage_level2_Fragments.parquet

Outputs (Analysis_Outputs/Lineage_Innervation-Focusness/):
    Lineage_Innervation_TBars_{level}.parquet
    Lineage_Innervation_PSDs_{level}.parquet
    Lineage_Innervation_Dendrites-PSDs_{level}.parquet
    Lineage_Focusness_{level}.parquet
    Lineage_Focusness_{level}_Dendrites.parquet
"""
import numpy as np
import pandas as pd
import os

from scipy.cluster.hierarchy import leaves_list, linkage

#------------------------------------------------------------------------------
FOLDER = "Analysis_Outputs/Lineage_Innervation-Focusness"
HOME   = os.environ["PROJECTS_HOME"]

EXCLUDE_PATTERNS = ("S-", "Unbound", "BolwigsAxon")

#------------------------------------------------------------------------------
def load_synapse_by_lineage(level, fragments = False):
    """
    Read the Synapse-by-Lineage parquet from L1_Neuron-Typing.

    """
    suffix = "_Fragments" if fragments else ""
    path   = f"{HOME}/L1_Neuron-Typing/Analysis_Outputs/Synapses/Synapse-by-Lineage_{level}{suffix}.parquet"
    return pd.read_parquet(path)

#------------------------------------------------------------------------------
def filter_lineages(df):
    """
    Drop sensory, unbound and Bolwig's-axon entries.

    """
    mask = ~df["Lineage"].str.contains("|".join(EXCLUDE_PATTERNS))
    return df.loc[mask]

#------------------------------------------------------------------------------
def add_relative_compartment(df):
    """
    Collapse left/right hemispheres into Ipsi/Contra.

    Compartment names follow ``<name>_<Left|Right>`` (e.g. ``AL_Left``,
    ``MB-CA_Right``, ``LALcd_Right``); the bare name comes from the rsplit.

    """
    df = df.copy()

    bare = df["Compartment"].str.rsplit("_", n = 1).str[0]
    side = df["Compartment"].str.rsplit("_", n = 1).str[1]

    df["CompRelative"] = np.where(df["Hemisphere"] == side, "Ipsi", "Contra")
    df["Comp"]         = bare + "_" + df["CompRelative"]
    return df

#------------------------------------------------------------------------------
def shannon_entropy(row):
    """
    Shannon entropy of a normalized row (zeros ignored).

    """
    row = np.asarray(row, dtype = float)
    if not np.isclose(row.sum(), 1.0):
        row = row / row.sum()
    nz = row[row > 0]
    return -(nz * np.log(nz)).sum()

def simpson_index(row):
    """
    Herfindahl-Hirschman / Simpson index of a normalized row.

    """
    row = np.asarray(row, dtype = float)
    if not np.isclose(row.sum(), 1.0):
        row = row / row.sum()
    return float(np.sum(row ** 2))

#------------------------------------------------------------------------------
def build_matrix(df, lineages):
    """
    Pivot to (Lineage x Comp), average across hemispheres and row-normalize.

    Lineages without any synapses for this connector are kept as all-NaN
    rows after row-normalization (0 / 0). DPLc2 is appended explicitly so
    it remains visible in the saved matrices even when it has no data.

    """
    avg = (df.groupby(["Lineage", "Comp"])["Count"]
             .mean()
             .reset_index())

    matrix = avg.pivot(index = "Lineage", columns = "Comp", values = "Count")
    matrix = matrix.reindex(index = lineages).fillna(0)

    # Row-normalize; all-zero rows become all-NaN
    matrix = matrix.div(matrix.sum(axis = 1), axis = 0)
    return matrix

#------------------------------------------------------------------------------
def ward_reorder(*matrices):
    """
    Ward-cluster the row-fusion of the given matrices and reapply the
    leaf order (reversed, matching the notebooks) to each one.

    """
    fused = pd.concat(matrices, axis = 1, join = "inner").fillna(0)
    order = leaves_list(linkage(fused, method = "ward", metric = "euclidean"))[::-1]

    return [m.loc[fused.index[order]] for m in matrices]

#------------------------------------------------------------------------------
def focusness_metrics(matrix, prefix):
    """
    Apply Shannon / Simpson and derived metrics row-wise. Lineages whose
    rows are entirely NaN (no synapses for this connector) are dropped.

    """
    matrix = matrix.dropna(how = "all")
    out = pd.DataFrame({
        f"{prefix}_Simpson" : matrix.apply(simpson_index, axis = 1),
        f"{prefix}_Entropy" : matrix.apply(shannon_entropy, axis = 1),
    })
    out[f"{prefix}_Entropy_Neff"] = np.exp(out[f"{prefix}_Entropy"])
    out[f"{prefix}_Simpson_Neff"] = 1 / out[f"{prefix}_Simpson"]
    out[f"{prefix}_Hierarchy"]    = out[f"{prefix}_Entropy_Neff"] / out[f"{prefix}_Simpson_Neff"]
    return out

#------------------------------------------------------------------------------
def whole_pipeline(level):
    """
    Whole-neuron T-bar + PSD matrices and focusness for a given level.

    """
    df = load_synapse_by_lineage(level, fragments = False)
    df = filter_lineages(df)
    df = add_relative_compartment(df)

    # Keep DPLc2 visible as an empty row even though it has no synapses here
    lineages = np.append(np.unique(df["Lineage"]), "DPLc2")

    tbar_matrix = build_matrix(df.loc[df["Connector"] == "TBars"], lineages)
    psd_matrix  = build_matrix(df.loc[df["Connector"] == "PSDs"],  lineages)

    # Fused Ward clustering — same row order on both matrices
    tbar_matrix, psd_matrix = ward_reorder(tbar_matrix, psd_matrix)

    tbar_metrics = focusness_metrics(tbar_matrix, "TBars")
    psd_metrics  = focusness_metrics(psd_matrix,  "PSDs")

    # Outer join keeps lineages that have only one connector type
    # (e.g. TRvm — PSDs but no TBars) with NaN on the missing side
    focus = pd.concat((tbar_metrics, psd_metrics), axis = 1, join = "outer")
    focus["Breadth_Ratio"]   = focus["TBars_Entropy_Neff"] / focus["PSDs_Entropy_Neff"]
    focus["Hierarchy_Ratio"] = focus["TBars_Hierarchy"]    / focus["PSDs_Hierarchy"]

    return tbar_matrix, psd_matrix, focus.reset_index()

#------------------------------------------------------------------------------
def dendrite_pipeline(level):
    """
    Dendrite-only PSD matrix and focusness for a given level.

    """
    df = load_synapse_by_lineage(level, fragments = True)
    df = df.loc[(df["Piece"] == "dendrite") & (df["Connector"] == "PSDs")]
    df = filter_lineages(df)
    df = add_relative_compartment(df)

    lineages = np.append(np.unique(df["Lineage"]), "DPLc2")

    psd_matrix, = ward_reorder(build_matrix(df, lineages))

    return psd_matrix, focusness_metrics(psd_matrix, "PSDs").reset_index()

#------------------------------------------------------------------------------
if __name__ == "__main__":
    # Go through each level
    for level in ["level1", "level2"]:
        print(f"Level : {level}\n")

        # Whole-neuron matrices and focusness
        tbar_matrix, psd_matrix, whole_focus = whole_pipeline(level)
        tbar_matrix.reset_index().to_parquet(f"{FOLDER}/Lineage_Innervation_TBars_{level}.parquet", index = False, engine = "pyarrow")
        psd_matrix.reset_index().to_parquet(f"{FOLDER}/Lineage_Innervation_PSDs_{level}.parquet", index = False, engine = "pyarrow")
        whole_focus.to_parquet(f"{FOLDER}/Lineage_Focusness_{level}.parquet", index = False, engine = "pyarrow")

        # Dendrite-only PSD matrix and focusness
        dendrite_matrix, dendrite_focus = dendrite_pipeline(level)
        dendrite_matrix.reset_index().to_parquet(f"{FOLDER}/Lineage_Innervation_Dendrites-PSDs_{level}.parquet", index = False, engine = "pyarrow")
        dendrite_focus.to_parquet(f"{FOLDER}/Lineage_Focusness_{level}_Dendrites.parquet", index = False, engine = "pyarrow")

    print("हो गया दोस्तों!")
