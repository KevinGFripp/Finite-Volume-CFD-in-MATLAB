function L = PressureLaplacian(Mesh)
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

Lx = Laplacian_Nonuniform_Neumann(Mesh.Nx,Mesh.P_dx_centre,Mesh.dx);
Ly = Laplacian_Nonuniform_Neumann(Mesh.Ny,Mesh.P_dy_centre,Mesh.dy);

Iy = speye(Mesh.Ny,Mesh.Ny);
Ix = speye(Mesh.Nx,Mesh.Nx);

L = kron(Ly,Ix) + kron(Iy,Lx);

L(1,:) = 0;
L(1,1) = 1;

end
