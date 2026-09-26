# L1 Compartment Connectivity

Neuropil compartment wiring, synaptic density, and lineage innervation in the
first-instar *Drosophila* larval (L1) connectome. Characterizes how the ~60
brain compartments connect to each other and how individual neuron lineages
distribute their synapses across the neuropil.

## Overview

The pipeline aggregates per-neuron synapse tables (from `L1_Neuron-Typing`) and
compartment volumes (from `L1_Compartment-Atlas`) into compartment-level
connectivity matrices and lineage innervation profiles. Core outputs:

- **Compartment Connectivity** — directed, weighted wiring diagram between all
  Level 1 and Level 2 compartments, normalized by presynaptic (PSD) and
  postsynaptic (T-bar) statistics
- **Synaptic Density** — T-bar and PSD counts normalized by compartment volume,
  revealing hot spots of input and output connectivity
- **Lineage Innervation** — (Lineage × Compartment) innervation matrices for
  central brain and VNC lineages, separated by axonal outputs (T-bars),
  dendritic inputs (PSDs), and dendrite-only views
- **Lineage Focusness** — Shannon entropy, Simpson/HHI index, and breadth /
  hierarchy ratios quantifying how narrowly each lineage targets the neuropil
- **Sensory Innervation** — specialized innervation analysis for S-* and
  BolwigsAxon lineages at Level 2 resolution

## Pipeline

Scripts are intended to run sequentially; each builds on prior outputs.

| # | Script | Purpose |
|---|--------|---------|
| 1 | `compartment_connectivity.py` | Aggregate per-neuron synapse tables into pairwise (source → target) compartment connectivity weights at Level 1 and Level 2. See [Connectivity Weight](#connectivity-weight) below. |
| 2 | `compartment_density.py` | Compute T-bar and PSD synaptic densities per compartment, normalized by compartment volume from `L1_Compartment-Atlas` STL meshes. |
| 3 | `lineage_focusness.py` | Build (Lineage × Compartment) innervation matrices for non-sensory lineages at Level 1 and Level 2 (T-bars, PSDs, dendrite-PSDs); compute focusness metrics. See [Focusness Metrics](#focusness-metrics) below. Dendrite view uses PSD-only; Ward clustering applied for row reordering. |
| 4 | `sensory_innervation.py` | Specialized (Lineage × Compartment) analysis for S-* and BolwigsAxon lineages at Level 2; row-normalized with Ward clustering. |

Shared utilities:
- `_compartment_palette.py` — canonical compartment orderings and color palettes
  (`COMP_HEX`, `COMP_ORDER`) used across all figures and notebooks.
- `viz/heatmap.py` — reusable `plot_connectivity_heatmap()` for Ward-clusterable,
  normalizable connectivity heatmaps.

## Connectivity Weight

For a source compartment $s$ and target compartment $t$, the weight is defined
over the set of neurons $\mathcal{S}$ shared between the two compartments (i.e.
neurons with at least one PSD in $s$ and at least one T-bar fragment in $t$):

$$w(s \to t) = \left(\sum_{i \in \mathcal{S}} \text{PSDs}_i(s)\right) \times \frac{\displaystyle\sum_{i \in \mathcal{S}} \text{TBars}_i(t)}{\displaystyle\sum_{j \in t} \text{TBars}_j(t)}$$

The first term is the total postsynaptic count in the source compartment
contributed by shared neurons — a measure of how much input they receive there.
The second term is the fraction of total output in the target compartment that
those same neurons account for, normalizing by the **total number of T-bars
across all neurons in the target compartment**. Together the weight captures
how strongly neurons that are postsynaptic in $s$ collectively drive output
in $t$, normalized for the overall synaptic load of $t$.

## Focusness Metrics

Each row of a (Lineage × Compartment) innervation matrix is row-normalized to
sum to 1, giving a probability distribution $\mathbf{p}$ over compartments.
Focusness metrics are then computed from $\mathbf{p}$ per lineage:

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| Shannon Entropy | $H = -\sum_i p_i \ln p_i$ | Overall spread; 0 = one compartment, $\ln N$ = uniform |
| $N_\text{eff}^H$ | $\exp(H)$ | Effective number of compartments under entropy |
| Simpson Index | $D = \sum_i p_i^2$ | Dominance; 1 = one compartment, $1/N$ = uniform |
| $N_\text{eff}^D$ | $1/D$ | Effective number of compartments under Simpson |
| Hierarchy | $N_\text{eff}^H \;/\; N_\text{eff}^D$ | Shape of the distribution: 1 = uniform, $> 1$ = one dominant compartment with a tail |
| Breadth Ratio | $N_\text{eff}^H(\text{TBars}) \;/\; N_\text{eff}^H(\text{PSDs})$ | Whether axonal outputs span more compartments than dendritic inputs |
| Hierarchy Ratio | $\text{Hierarchy}(\text{TBars}) \;/\; \text{Hierarchy}(\text{PSDs})$ | Whether outputs are more hierarchically concentrated than inputs |

$N_\text{eff}^H$ and $N_\text{eff}^D$ both measure the effective number of
compartments a lineage innervates, but weight the tail differently: $N_\text{eff}^H$
is sensitive to rare compartments, while $N_\text{eff}^D$ is dominated by the
strongest ones. Their ratio (Hierarchy) therefore captures distribution shape
independently of breadth — a lineage with one strong target and many weak ones
has high Hierarchy, while a lineage spread evenly across a few compartments has
Hierarchy near 1.

## Figures

Each notebook generates a set of main and supplementary panels saved to `Figures/`.

| Notebook | Output |
|----------|--------|
| `Figure_02_Compartment-Connectivity.ipynb` | Synapse cloud projection, LAL anatomy, Level 1 connectivity heatmap, Level 2 brain / VNC connectivity (5 supplementary panels) |
| `Figure_02_Compartment-Synaptic-Density.ipynb` | Volume-normalized synaptic density per compartment |
| `Figure_03_Lineage-Innervation-Level1.ipynb` | Lineage × compartment innervation at Level 1; threshold sensitivity supplementary |
| `Figure_03_Lineage-Innervation-Level2.ipynb` | Lineage × compartment innervation at Level 2; T-bars, PSDs, and VNC segments |
| `Figure_03_Sensory-Innervation.ipynb` | Sensory lineage innervation heatmaps (T-bars + PSDs) |

## Directory Structure

```
L1_Compartment-Connectivity/
  Scripts/                            # Analysis scripts (+ viz/ subfolder)
  Notebooks/                          # Figure generation notebooks
  Data_Helpers/                       # Reference orderings, palettes, and morphologies
    Compartment_Colors_Level1.pkl
    Compartment_Order_*.npy
    Lineages_*.npy
    LAL-Neurons/                      # SWC morphology files
    Tracts/                           # STL mesh files
  Analysis_Outputs/
    Connectivity/                     # Compartment connectivity matrices (Level 1 & 2)
    Density/                          # Synaptic density tables (Level 1 & 2)
    Lineage_Innervation-Focusness/    # Innervation matrices + focusness metrics
    Sensory_Innervation/              # Sensory lineage innervation matrices
  Figures/
    Figure_02/                        # Compartment connectivity & density figures
    Figure_03/                        # Lineage innervation figures
```
