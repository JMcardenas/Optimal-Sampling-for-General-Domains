# Optimal-Sampling-
Solve weighted-LSP approximation of a function f on an irregular domain D, from samples points under a optimal measures.
* Function: Opt_sam_v1.m 
  This function sampling points in any iteration with respect the dimension and polynomial degree.
  Basic sctructure: 
  - Set up: Define the parameters. 
  - Generate K points under uniform measure: Generate the grid. 
  - Pre-processing:
    - Under hyperbolic cross index this generate n_max and N_max (polynomial degree and number of bases). 
    - Reorder index set: The first index set correspond to n_1 and the second index until n_2 and continues.
    - Generate B matrix: This matrix contain the Legendre polynomials evaluates on the grid. 
    - Compute measure mu.
    - Evaluate f on the grid.
  - Compute W-LS approximation:
    
          
