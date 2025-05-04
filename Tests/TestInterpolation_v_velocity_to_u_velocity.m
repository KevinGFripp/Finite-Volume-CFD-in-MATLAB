function TestInterpolation_v_velocity_to_u_velocity()

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
GrowthRatesx = 0.002;
GrowthRatesy = 0.002;

Mesh = MakeRectilinearMesh(Nx,Ny,Width,Height,GrowthRatesx,GrowthRatesy,[]);

% u velocity data
u_velocity = ones((Nx+1)*Ny,1).*linspace(1,(Nx+1)*Ny,(Nx+1)*Ny).';


% v velocity data
v_velocity = ones(Nx*(Ny+1),1);

v_x_Grid = repmat(Mesh.v_centre_x,Ny+1,1);
v_y_Grid = repmat(Mesh.v_centre_y,Nx,1);

u_x_Grid =repmat(Mesh.u_centre_x,Ny,1);
u_y_Grid = repmat(Mesh.u_centre_y,Nx+1,1);

u_Mesh_Interpolant = scatteredInterpolant(u_x_Grid,u_y_Grid,ones(Nx*(Ny+1),1));
v_Mesh_Interpolant = scatteredInterpolant(v_x_Grid,v_y_Grid,ones((Nx+1)*Ny,1));

u_Mesh_Interpolant.Values = u_velocity;
v_Mesh_Interpolant.Values = v_velocity;

u_interp = u_Mesh_Interpolant(v_x_Grid,v_y_Grid);
v_interp = v_Mesh_Interpolant(u_x_Grid,u_y_Grid);


imagesc(reshape(u_interp,Nx,Ny+1))



end