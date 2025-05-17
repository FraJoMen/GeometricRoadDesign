# Geometric Road Design

This project explores the mathematical and computational aspects of road geometric design, focusing on:

- **Planimetric curve transitions** (horizontal alignments)  
- **Altimetric curve transitions** (vertical alignments)  
- **Comprehensive road design verifications**  

Jupyter Notebooks provide a didactic framework for exploring theoretical concepts and practical applications in road alignment analysis and design.

The /notebooks folder contains Jupyter Notebooks developed for road geometric design.

## Notebooks

This project is a **work in progress**, and the `/notebooks` folder may not yet contain all the planned files. The current and expected Jupyter Notebooks include:  

- **Planimetric Curve Transitions**: Study of circular arcs, clothoids, and their representation in a Cartesian plane.
- **Altimetric Curve Transitions** *(in progress)*: Analysis of vertical curves (concave and convex transitions).  
- **Advanced Curve Transitions** *(planned)*: More complex configurations involving consecutive curves in the same or opposite direction, including flexure clothoids and minimum separation conditions.
- **Road Design Verification** *(planned)*: Computational checks for road geometry compliance.  
  

## Installation

This project is based on Jupyter Notebooks and requires no formal installation.  
We recommend using the [Anaconda distribution](https://www.anaconda.com/) to ensure access to the necessary scientific Python libraries such as `numpy`, `plotly`, and `matplotlib`.

To use the project locally:

1. Clone or download this repository.
2. Open the desired notebook from the `/notebooks` folder in your local Jupyter environment.

All notebooks are designed to run independently and provide a self-contained environment for exploring geometric road design concepts.

## Python Modules

The core functionalities for geometric road alignment computations are implemented in the `functions.py` module, which defines reusable classes and tools for use across multiple notebooks.

Currently implemented:

- `CircularTransition`: handles geometric transitions between straight segments using circular arcs and clothoids. Includes:
  - Geometry construction (intersection angles, circle center, tangent points)
  - Clothoid generation and curvature correction
  - Visualizations using Plotly
  - Convergence analysis for clothoid series expansions
  

## Contributing

Contributions are welcome! Feel free to fork the repository, open issues, or submit pull requests.  

## License

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.  
