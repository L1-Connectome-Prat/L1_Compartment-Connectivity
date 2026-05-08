#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-compartment T-bar and PSD synaptic densities. Synapse counts come from
the per-neuron synapse tables produced by L1_Neuron-Typing; compartment
volumes come from the STL meshes in L1_Compartment-Atlas.

Inputs (read-only):
    $PROJECTS_HOME/L1_Neuron-Typing/Analysis_Outputs/Synapses/
        Synapse-by-Neuron_level1.parquet
        Synapse-by-Neuron_level2.parquet
    $PROJECTS_HOME/L1_Compartment-Atlas/
        Level1.npy
        Level2.npy

Outputs (Analysis_Outputs/Density/):
    Compartment_Density_level1.parquet
    Compartment_Density_level2.parquet
"""
import numpy as np
import pandas as pd
import os

import connectome_analysis_claude as ca

#------------------------------------------------------------------------------
SCALE  = 1e9
FOLDER = "Analysis_Outputs/Density"
HOME   = os.environ["PROJECTS_HOME"]
COMPARTMENTS_HOME = f"{HOME}/L1_Compartment-Atlas"

#------------------------------------------------------------------------------
def load_synapse_by_neuron(level):
    """
    Read the Synapse-by-Neuron parquet from L1_Neuron-Typing.

    """
    path = f"{HOME}/L1_Neuron-Typing/Analysis_Outputs/Synapses/Synapse-by-Neuron_{level}.parquet"
    return pd.read_parquet(path)

#------------------------------------------------------------------------------
def load_compartment_volumes(level):
    """
    Load compartment STL meshes from L1_Compartment-Atlas and return a
    {name : volume} dictionary.

    """
    compartment_file  = f"{COMPARTMENTS_HOME}/{level.capitalize()}.npy"
    compartment_paths = np.load(compartment_file, allow_pickle = True)

    volumes = {}
    for comp in compartment_paths:
        stl       = ca.stl_to_navis(f"{COMPARTMENTS_HOME}/{comp}")
        name      = comp.split("/")[-1].replace(".stl", "")
        volumes[name] = stl.volume

    return volumes

#------------------------------------------------------------------------------
def compute_density(neuron_df, volumes):
    """
    Aggregate synapse counts per compartment and divide by mesh volume.

    """
    summed = (neuron_df
              .groupby(["Compartment", "Connector"])["Count"]
              .sum()
              .unstack(fill_value = 0)
              .rename(columns = {"TBars" : "Num_TBars", "PSDs" : "Num_PSDs"})
              .reset_index())

    summed["Volume"] = summed["Compartment"].map(volumes)
    summed = summed.dropna(subset = ["Volume"])

    summed["Vol_Norm"]      = summed["Volume"] / summed["Volume"].max()
    summed["TBar_Vol_Norm"] = (summed["Num_TBars"] / summed["Volume"]) * SCALE
    summed["PSD_Vol_Norm"]  = (summed["Num_PSDs"]  / summed["Volume"]) * SCALE

    return summed

#------------------------------------------------------------------------------
if __name__ == "__main__":
    # Go through each level
    for level in ["level1", "level2"]:
        print(f"Level : {level}\n")

        # Load synapse counts and compartment volumes
        neuron_df = load_synapse_by_neuron(level)
        volumes   = load_compartment_volumes(level)

        # Compute density
        density_df = compute_density(neuron_df, volumes)
        density_df.to_parquet(f"{FOLDER}/Compartment_Density_{level}.parquet", index = False, engine = "pyarrow")

    print("हो गया दोस्तों!")
