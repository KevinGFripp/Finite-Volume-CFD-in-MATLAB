function TestInterpolation_u_Velocity_To_vgrid_PressureFaces()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

Width = 1.0;
Height = 1.0;
Nx = 32;
Ny = 32;
GrowthRatesx = 2;
GrowthRatesy = 2;

Mesh = MakeRectilinearMesh(Nx,Ny,Width,Height,GrowthRatesx,GrowthRatesy,[]);

Ru = Interpolate_u_velocity_to_face(Mesh);

u_velocity = ones((Nx+1)*Ny,1).*linspace(0,1,(Nx+1)*Ny).';

figure(1)
imagesc(reshape(Ru*u_velocity,Nx,Ny+1))

figure(2)
imagesc(full(Ru))

ExpectedSum = length(nonzeros(Ru))/2 +3/2 * Nx
ActualSum = full(sum(Ru,"all"))


end

function k = index(i,j,Nx)

k = i +(j-1)*Nx;

end