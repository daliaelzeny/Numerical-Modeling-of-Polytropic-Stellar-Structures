# Numerical-Modeling-of-Polytropic-Stellar-Structures
Numerical solution of the Lane-Emden equation for polytropic indices n = 1.5, 3, 4 in Python, scaled to physical units and compared with MESA stellar models and Gaia data.
Course project for **PHYS 331: Stellar Structure and Evolution**, Zewail University of
Science and Technology, May 2026.

**Authors:** Dalia Zakaria, Mirna Ahmed

## Overview

A polytropic star follows the equation of state P = Kρ^(1 + 1/n). Combining it with
hydrostatic equilibrium and mass conservation gives the Lane-Emden equation:

    (1/ξ²) d/dξ ( ξ² dθ/dξ ) = −θⁿ,    θ(0) = 1,  θ'(0) = 0

where ρ = ρc θⁿ and r = αξ.

This project:

1. Derives the Lane-Emden equation analytically from hydrostatic equilibrium and mass conservation
2. Integrates it numerically in Python (Euler method) for n = 1.5, 3 and 4
3. Scales the dimensionless solutions to physical units (M = M☉, R = R☉) to obtain
   α, ρc, Pc, Tc and the density, pressure and temperature profiles
4. Compares the normalized density profiles with a 1 M☉ ZAMS model from MESA, using mean squared error (MSE)
5. Compares the profiles qualitatively with stellar populations seen in Gaia
   (low-mass convective stars, white dwarfs, solar-type stars)

## Results

- **n = 3** matches the MESA solar-type profile best (MSE ≈ 0.0007), consistent with
  the Eddington standard model for a star dominated by radiative energy transport.
- **n = 1.5** (fully convective) decreases too slowly compared with the MESA profile.
- **n = 4** is too centrally condensed and leaves an almost empty envelope, so it is
  not a realistic model of a stable main-sequence star.

## Repository structure

    ├── code/
    │   └── lane_emden.py       # solver, scaling, MESA comparison, plots
    ├── data/
    │   └── zams_1.0.mod        # MESA model (see Data below)
    ├── plots/                  # θ(ξ), density, pressure, temperature, MESA comparison
    ├── report.pdf              # full project report
    └── README.md

## Usage

    pip install numpy scipy matplotlib
    python code/lane_emden.py

The script prints ξ₁, dθ/dξ at the surface, α, ρc, Pc, Tc and the MSE for each n,
then shows the plots.

## Data

The MESA comparison uses the 1.0 M☉ ZAMS model file `zams_1.0.mod` from:

Li, B. et al. (2025). *MESA files for the binary estimations of M2C*. Zenodo.
https://doi.org/10.5281/zenodo.15778913

## Method notes

- Integration starts at ξ = 10⁻⁶ with step 10⁻⁵ (Euler method) and stops when θ < 0;
  the surface ξ₁ is found by linear interpolation.
- Both the Lane-Emden and MESA density profiles are normalized (r/R and ρ/ρc) and
  interpolated onto the MESA radial grid before computing the MSE.
- Limitations: the polytropic model assumes spherical symmetry, hydrostatic
  equilibrium, a single polytropic equation of state and constant composition. MESA
  includes nuclear burning, realistic opacities and equations of state, and
  composition changes, so some deviation is expected.
- The Gaia comparison is qualitative: it relates the model profiles to the stellar
  populations seen in Gaia HR diagrams, not a fit to Gaia data.

## Contributions

- **Dalia Zakaria:** numerical solution and Python code, MESA comparison, abstract, conclusion
- **Mirna Ahmed:** analytical derivation, Gaia comparison, introduction, references

## References

1. Prialnik, D. (2009). *An Introduction to the Theory of Stellar Structure and
   Evolution*. Cambridge University Press.
2. Paxton, B. et al. (2011). MESA: Modules for Experiments in Stellar Astrophysics.
   *ApJS*, 192(1), 3.
3. Li, B. et al. (2025). MESA files for the binary estimations of M2C. Zenodo.
   https://doi.org/10.5281/zenodo.15778913
4. Gaia Collaboration (2018). *A&A*, 616, A1.
