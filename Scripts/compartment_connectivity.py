#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute pairwise source → target compartment connectivity weights from the
per-neuron synapse tables produced by L1_Neuron-Typing.

Weight definition (per compartment pair):

    w(s -> t) = (Σ PSDs of shared neurons in source)
              × (Σ T-bars of shared neurons in target)
              / (Σ T-bars of all fragments in target)

Source neurons are required to have more than ``PSD_THRESH`` PSDs in their
source compartment (matches the original ``>= 1`` filter).

Inputs (read-only):
    $PROJECTS_HOME/L1_Neuron-Typing/Analysis_Outputs/Synapses/
        Synapse-by-Neuron_level1.parquet
        Synapse-by-Neuron_level2.parquet

Outputs (Analysis_Outputs/Connectivity/):
    Compartment_Connectivity_level1.parquet
    Compartment_Connectivity_level2.parquet
"""
import pandas as pd
import os

#------------------------------------------------------------------------------
PSD_THRESH = 1
FOLDER = "Analysis_Outputs/Connectivity"
HOME   = os.environ["PROJECTS_HOME"]
#------------------------------------------------------------------------------
def load_synapse_by_neuron(level):
    """
    Read the Synapse-by-Neuron parquet from L1_Neuron-Typing.
    
    """
    path = f"{HOME}/L1_Neuron-Typing/Analysis_Outputs/Synapses/Synapse-by-Neuron_{level}.parquet"
    return pd.read_parquet(path)

#------------------------------------------------------------------------------
def neuron_table_to_pivot(neuron_df):
    """
    Long table → DataFrame with columns [Compartment, bodyId, TBars, PSDs].

    """
    pivot = (neuron_df
             .pivot_table(index   = ["Compartment", "bodyId"],
                          columns = "Connector",
                          values  = "Count",
                          aggfunc = "sum")
             .fillna(0)
             .reset_index())

    for col in ("TBars", "PSDs"):
        if col not in pivot.columns:
            pivot[col] = 0

    return pivot[["Compartment", "bodyId", "TBars", "PSDs"]]


#------------------------------------------------------------------------------
def pairwise_compartment_connectivity(pivot, psd_thresh = PSD_THRESH):
    """Pairwise compartment connectivity weights from a per-(comp, neuron) pivot."""
    comps  = pivot["Compartment"].unique()
    by_cmp = {c : g.set_index("bodyId") for c, g in pivot.groupby("Compartment")}

    rows = []

    for tgt in comps:
        tgt_frags   = by_cmp[tgt]
        total_tbars = tgt_frags["TBars"].sum()

        if total_tbars == 0:
            continue

        for src in comps:
            src_frags = by_cmp[src]
            src_frags = src_frags[src_frags["PSDs"] >= psd_thresh]

            shared = src_frags.index.intersection(tgt_frags.index)
            if len(shared) == 0:
                continue

            source_psds  = src_frags.loc[shared, "PSDs"].sum()
            target_tbars = tgt_frags.loc[shared, "TBars"].sum()

            rows.append({
                "Source_Compartment" : src,
                "Target_Compartment" : tgt,
                "Total_TBars"        : total_tbars,
                "Target_TBars"       : target_tbars,
                "Source_PSDs"        : source_psds,
                "Weight"             : source_psds * (target_tbars / total_tbars)})

    return pd.DataFrame(rows)


#------------------------------------------------------------------------------
if __name__ == "__main__":
    # Go through each level
    for level in ["level1", "level2"]:
        print(f"Level : {level}\n")

        # Load synapse by neuron data
        neuron_df = load_synapse_by_neuron(level)
        # Pivot the data to get the desired format
        pivot     = neuron_table_to_pivot(neuron_df)

        weight_df = pairwise_compartment_connectivity(pivot)
        weight_df.to_parquet(f"{FOLDER}/Compartment_Connectivity_{level}.parquet", index  = False, engine = "pyarrow")

    print("हो गया दोस्तों!")
