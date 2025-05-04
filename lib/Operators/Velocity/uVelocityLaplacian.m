function [Lu,bu_bc] = uVelocityLaplacian(Mesh,Boundaries)
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%


%% velocity defined on left and right boundaries
 Lux = Laplacian_Nonuniform(Mesh.Nx+1,Mesh.u_dx_centre,Mesh.u_dx);


%% velocity boundary conditions top and bottom
[Luy,by_bc] = Laplacian_Nonuniform_Dirichlet(Mesh.Ny, ...
                                     Mesh.u_dy_centre, ...
                                     Mesh.u_dy, ...
                                     Boundaries.Bottom_u, ...
                                     Boundaries.Top_u);

Iy = speye(Mesh.Ny,Mesh.Ny);
Ix = speye(Mesh.Nx+1,Mesh.Nx+1);

Lu = kron(Luy,Ix) + kron(Iy,Lux);

bu_bc = zeros((Mesh.Nx+1)*Mesh.Ny,1);

BottomWall = index(1:Mesh.Nx+1,1,Mesh.Nx+1);
TopWall = index(1:Mesh.Nx+1,Mesh.Ny,Mesh.Nx+1);

for n= 1:length(BottomWall)
   bu_bc(BottomWall(n)) = by_bc(1);
   bu_bc(TopWall(n)) = by_bc(Mesh.Ny);
end





end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end