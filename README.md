# Finite-Volume CFD in MATLAB
A 2D time-dependent finite-volume MATLAB solver for the Navier-Stokes equations on non-uniform cartesian meshes with staggered grids.

# Implementation
- Uniform and non-uniform Cartesian 2D staggered grids.
- Neumann and Dirichlet boundary conditions on grid boundaries or rectangular shapes.
- All differential operators are 2nd order accurate.
- Incompressibility is enforced by a pressure correction Poisson solve.
- Explicit and Implicit Runge-Kutta methods of different orders can be used to perform time-integration on the Navier-Stokes equations with adaptive step size.

# Validation : Lid Driven Cavity Re = 1000

<img width="892" height="462" alt="LidDrivenCavityResult" src="https://github.com/user-attachments/assets/af310701-c97e-4a34-acb9-6e698c55b35a" />


# Example : Flow Past a Square Cylinder Re = 120

https://github.com/user-attachments/assets/0b842add-72db-4b05-b62a-b2a9df807f08


# Example : Lid Driven Cavity Re = 10000

https://github.com/user-attachments/assets/6e195455-a07c-4b67-a426-12304e55dab2



