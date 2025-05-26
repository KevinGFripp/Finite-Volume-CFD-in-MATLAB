A 2D time-dependent finite-difference MATLAB solver for the Navier-Stokes equations on non-uniform cartesian meshes with staggered grids.

Two example simulations are included:
Lid driven cavity at Re = 1000
Flow past a square cylinder at Re = 120


The spatial derivatives are computed with 2nd order accuracy.

Time integration is performed with adaptive-stepping of either explicit Runge-Kutta (3rd or 5th order) or Implicit ESDIRK methods (3rd or 5th order).
