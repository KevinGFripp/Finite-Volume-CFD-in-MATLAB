function PoissonSolution_RegularGridReference()
%% include lib
basefolder = ...
'C:\Users\kevin\Documents\MANNGA\Data\April 25\2D_StaggeredGrid_NonUniformMesh_FiniteVolume\';
addpath(genpath(basefolder));
%%

% Grid
Nx = 80;
Ny = 80;
W = 1;
H = 1;

dx = W/Nx;
dy = H/Ny;

[Lxy,b] = LaplacianMatrix_RegularGrid(Nx,Ny,dx,dy,-1.0);

t = Lxy\b;

imagesc(reshape(t,Nx,Ny))

end
