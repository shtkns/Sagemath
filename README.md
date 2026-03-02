# SageMath: Magnitude-Path Spectral Sequence (MPSS)

A SageMath implementation for computing **Magnitude-Path Spectral Sequences** of graphs, supporting up to the 2nd sheet (Bigraded Path Homology).

## Sample Code

Use this snippet to quickly verify the functionality and see how the spectral sequence sheets are generated.

```python
from sage.all import *
# Ensure mpss_tools.py is in your directory
from mpss_tools import MPSS 

# 1. Setup: Create a Cycle Graph (C_5)
g = graphs.CycleGraph(5)

# 2. Initialize MPSS up to the 2nd sheet
# 'side' sets the max range for both degree k and length l
E = MPSS(graph=g, side=6, rmax=2)

# 3. Output: Print the 1st sheet (Magnitude Homology)
print("--- Magnitude Homology (E1 Sheet) ---")
E.print_sheet(1)

# 4. Inspection: Retrieve specific generators
k, l, r = 2, 2, 1
gens = E.get_generators(k, l, r)
print(f"\nGenerators for (k={k}, l={l}) on sheet {r}:")
print(gens)
```

## Key Features

### Spectral Sequence Computation
The `MPSS` class automates the construction of chain complexes and transitions between sheets:
- **`print_sheet(r)`**: Generates a formatted rank table for the $r$-th sheet ($r=0, 1, 2$). Note that the table only displays the free ranks and **does not contain torsion information.**
- **`magnitude()`**: Computes the magnitude of the graph as a formal power series string, **truncated up to the degree specified by `side`.**

### Generator & Space Analysis
- **`get_generators(k, l, r)`**: Returns the actual cycles (paths/chains) representing the homology generators.
- **`get_divisors(k, l, r)`**: Retrieves the generators of the boundaries (images of the differential). This helps identify the relations in the homology group and analyze which classes become trivial.
- **`show_space_info(k, l, r)`**: Displays the algebraic structure of the specific $(k, l, r)$ space, including torsion information (e.g., $\mathbb{Z}^n \oplus \mathbb{Z}/m_1\mathbb{Z} \oplus \dots$). **(Note: This feature is currently a work in progress.)**

## Parameters of the `MPSS` Class

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `graph` | Sage Graph | *Required* | The SageMath graph object to be analyzed. |
| `side` | int | *Required* | Max boundary for bigrading indices $k$ and $\ell$. |
| `rmax` | int | `2` | Maximum sheet index ($0, 1, 2$) to be computed. |
| `Reduced` | bool | `True` | Whether to compute the reduced version of homology. |
| `Eulerian` | bool | `False` | If True, computation is restricted to Eulerian Magnitude Homology. |

## Installation
Download `mpss_tools.py` and place it in the same directory as your execution script.

```text
.
├── mpss_tools.py
└── your_script.sage (or .py)

```
