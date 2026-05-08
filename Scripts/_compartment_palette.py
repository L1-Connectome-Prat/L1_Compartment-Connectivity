"""
Compartment palette + ordering helpers used across the
compartment-connectivity figures.

Static data lives in ``Data_Helpers/`` at the repo root so the orderings,
colour map and node-bubble layout can be loaded independently (e.g. from
notebooks or sibling repos) without importing this module:

    Data_Helpers/
        Compartment_Order_Level1.npy
        Compartment_Order_Level2.npy
        Relative_Compartment_Order_Level1.npy
        Compartment_Colors_Level1.pkl
        NodeBubble.pkl

This module exposes the same data plus the high-level acronym palette
shared with L1_Network.
"""
import os
import pickle as pkl

import numpy as np

#------------------------------------------------------------------------------
DATA_DIR = os.path.join(os.environ["PROJECTS_HOME"],
                        "L1_Compartment-Connectivity",
                        "Data_Helpers")

#------------------------------------------------------------------------------
COMP_ORDER = ["AL", "TR", "LON",
              "MB-CA", "MB", "MBE",
              "IPA", "IPLM", "IPP",
              "SMP", "SMPal",
              "SLP",
              "LAL", "VMC", "LAL-VMC",
              "VLP",
              "SEZ-VNC",
              "Mixed"]


#------------------------------------------------------------------------------
def load_compartment_order(level = "level1", relative = False):
    """Return the canonical row/column ordering as a list of strings."""
    if relative:
        if level != "level1":
            raise ValueError("Relative ordering is only defined for level1.")
        name = "Relative_Compartment_Order_Level1.npy"
    else:
        name = {"level1" : "Compartment_Order_Level1.npy",
                "level2" : "Compartment_Order_Level2.npy"}[level]

    return np.load(os.path.join(DATA_DIR, name)).tolist()


def load_compartment_colors(level = "level1"):
    """Return the PK_*-keyed compartment colour dict."""
    name = {"level1" : "Compartment_Colors_Level1.pkl"}[level]
    with open(os.path.join(DATA_DIR, name), "rb") as f:
        return pkl.load(f)


def load_node_bubble():
    """Return the network bubble-layout dict (`position`, `volume`)."""
    with open(os.path.join(DATA_DIR, "NodeBubble.pkl"), "rb") as f:
        return pkl.load(f)
