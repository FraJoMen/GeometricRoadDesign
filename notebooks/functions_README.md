# functions.py – Core Module for Geometric Road Design

This module implements the computational core of the *GeometricRoadDesign* project.  
It provides reusable Python classes and methods for the geometric analysis of road alignments, focusing on planimetric transitions using circular arcs and clothoids.

---

## Class: `CircularTransition`

The main class `CircularTransition` models the transition between two straight segments using either a circular arc or a clothoid–arc–clothoid combination.

### 🔧 Geometric Construction

- Computes the **intersection angle** between incoming and outgoing directions.
- Determines the **center of the circular arc** and **tangent points**.
- Supports **center correction** based on clothoid geometry (O*).
- Identifies points T₁*, T₂*, C₁*, C₂* along the alignment.

### 📐 Clothoid Generation

- Supports **entry and exit clothoids** defined by parameters A₁ and A₂.
- Uses **series expansion** to generate clothoid points.
- Enforces **geometric continuity** at key locations (`T*`, `C*`).

### 📊 Visualization Tools

- High-quality **interactive plots** using Plotly.
- Displays road alignment segments, arc, clothoids, and key points.
- Configurable appearance and layout.

### 📈 Convergence Analysis

Includes utilities to analyze the convergence of series used in clothoid generation:

- `delta_R_convergence(A, R)`: evaluates offset between arc and clothoid.
- `Xm_convergence(A, R)`: midpoint displacement of the clothoid.
- `clothoid_convergence(A, R, τ)`: full convergence of (X(τ), Y(τ)).

---

## Example Usage

```python
from functions import CircularTransition

# Define three alignment points and a radius
P0 = [0, 0]
P1 = [50, 0]
P2 = [50, 50]
R = 80

# Initialize transition and add clothoids
trans = CircularTransition(P0, P1, P2, R)
trans.add_clothoids(A1=150, A2=150)

# Plot the result
trans.plotClothoid()
