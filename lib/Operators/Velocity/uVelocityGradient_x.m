function Du = uVelocityGradient_x(Mesh,Boundaries)
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%


%% velocity doesn't require boundaries
Dux = Gradient_Nonuniform(Mesh.Nx+1, ...
                                  Mesh.u_dx_centre, ...
                                  Mesh.u_dx);


Iy = speye(Mesh.Ny,Mesh.Ny);

Du = kron(Iy,Dux);


end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end