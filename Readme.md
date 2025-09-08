# Fully Vectorized Linear Operator Formulation for Inverse Scattering Problems

This repository contains the code and resources for the Master's Thesis: **"Fully Vectorized Linear Operator Formulation for Inverse Scattering Problems"**.

## Repository Structure
- **main.jl** and the **Modules** folder: Contain most of the code used for calculations.
- **Plots** folder: Contains additional plots not discussed in the thesis.
- **aluminum_plate_defect_Q_working.ipynb**: Jupyter notebook used to run Salvus simulations. Note: Requires a Salvus license. This is an older version from early in the thesis, missing some defects and tests.
- **dispersion.ipynb**: Contains parameters used with the `lambwaves` package for the discussion on dispersion curves. Documentation for the package can be found in the thesis references.

## Notes
- Currently, there is no data included for running the simulations. This repository serves as a reference for those interested in implementing the methods themselves.
- Excuse the messy internal structure. I am not a CS student and it shows.
- The calculation routines are mostly optimized.
