# pseudobands
Stochastic pseudobands (piggybacks on parabands/nscf)

Constructs stochastic linear combinations of Bloch states to speed up GW calculation 
(GW reduced from O(N^4) to ~O(N^2))

See /scripts/README for explanation/single use, and for convergence testing.

See /workflow for scripts that run many iterations of pseudobands at once. Useful for parallelizing convergence testing.

See workflow diagram for useage directions:

![alt text](./Screen%20Shot%202022-06-01%20at%201.45.51%20PM.png?raw=true)
