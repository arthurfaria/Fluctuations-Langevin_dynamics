# Brief description of the files:

1. stat_langevin.py
 	- Gaussian stochastic dynamics for a Brownian particle. Mean value and variance are calculated. Algorithm based on Ref.[1].
	
2. FT_langevin.py
	- Gaussian stochastic dynamics for a Brownian particle under a harmonic potential whose center of mass is displaced with constant velocity (see Ref.[2]). 
	Fluctuation Theorem for work is calculated. Algorithm based on Ref.[1].
	
## References:
[1] https://aip.scitation.org/doi/abs/10.1063/1.4802990
[2] https://journals.aps.org/pre/abstract/10.1103/PhysRevE.67.046102


# Quantum refrigerator

This is the repository for projects in collaboration with the NV centre group at the Institute of Physics III, University of Stuttgart, on theory and experimental implementations of quantum refrigerators.

## Repository content

1. Notes elaborating my personal understanding of the topics relevant in the research.
2. Overleaf git modules containing collaborative content, mainly of manuscript drafts.
3. Files for data analysis, plotting data, numerical and symbolic computation (such as raw `.csv` files, Mathematica notebooks `.nb`, packages `.m`, Jupyer notebooks `.ipynb`, etc).
4. File divisions for subprojects, such as those referring to the minimal HBAC, and those for the QET-HBAC.
5. Research references.
6. Files for presentations, posters, etc.

Disclaimer: not every file is being tracked by git and synchronized with remote repositories.

## Folder structure

The repository is divided into two projects, each concerning of the parts introduced in the previous section.

### Analysing the experimental data (Analysis)

1. `preliminary-run/`
