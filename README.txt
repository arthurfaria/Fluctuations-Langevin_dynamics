Brief files description:

1. stat_langevin.py
 	- stochastic dynamics for a Brownian particle. Mean value and variance are calculated. Algorithm based on Ref.[1]"
	
2. FT_langevin.py
	- stochastic dynamics for a Brownian particle under a harmonic potential whose center of mass is displaced with constant velocity (see Ref.[2]). 
  Fluctuation Theorem for work and heat are calculated. Algorithm based on Ref.[1]"
	
3. Generalized_FP.py
	- generalized Fokker-Planck (GenBM) equation (see 'gener_FP.pdf' file for further infos). Both the generalized semiclassical distriubtion (GenBM_rho) 
  and the distribution of a standard Brownian motion (BM_rho) are obatined considering a external harmonic potential. 
  Algorithm based on finite diference approach to compute derivatives.
	
References:
[1] https://aip.scitation.org/doi/abs/10.1063/1.4802990
[2] https://journals.aps.org/pre/abstract/10.1103/PhysRevE.67.046102
