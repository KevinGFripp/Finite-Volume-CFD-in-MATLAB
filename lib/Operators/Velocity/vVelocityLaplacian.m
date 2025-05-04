function [Lv,bv_bc] = vVelocityLaplacian(Mesh,Boundaries)
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%


%% velocity boundary conditions on left and right
[Lvx,bx_bc] = Laplacian_Nonuniform_Dirichlet(Mesh.Nx, ...
                                     Mesh.v_dx_centre, ...
                                     Mesh.v_dx,...
                                     Boundaries.Left_v, ...
                                     Boundaries.Right_v);

%% velocity defined on top and bottom
Lvy = Laplacian_Nonuniform(Mesh.Ny+1, ...
                           Mesh.v_dy_centre, ...
                           Mesh.v_dy);

Iy = speye(Mesh.Ny+1,Mesh.Ny+1);
Ix = speye(Mesh.Nx,Mesh.Nx);

Lv = kron(Lvy,Ix) + kron(Iy,Lvx);

bv_bc = zeros(Mesh.Nx*(Mesh.Ny+1),1);

LeftWall = index(1,1:Mesh.Ny+1,Mesh.Nx);
RightWall = index(Mesh.Nx,1:Mesh.Ny+1,Mesh.Nx);

for n= 1:length(LeftWall)
   bv_bc(LeftWall(n)) = bx_bc(1);
   bv_bc(RightWall(n)) = bx_bc(Mesh.Nx);
end


end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end