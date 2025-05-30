function Dv = vVelocityGradient_y(Mesh,Boundaries)
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%


%% velocity doesn't require boundary conditions
 Dvy = Gradient_Nonuniform(Mesh.Ny+1, ...
                           Mesh.v_dy_centre, ...
                           Mesh.v_dy);


Ix = speye(Mesh.Nx,Mesh.Nx);
Dv = kron(Dvy,Ix);


end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end