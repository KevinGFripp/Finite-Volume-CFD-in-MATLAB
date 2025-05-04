function PressureLaplacian_ShapeBoundaries()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

% Grid
Nx = 64;
Ny = 64;
Width = 1;
Height = 1;
% growth rate
wx = 0.5;
wy = 0.5;

Square = DefineBody(0.3,0.3,Width/Nx,Height/Ny,Nx,Ny,Width/Nx,Height/Ny);
Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height,[wx wx wx],[wy wy wy],Square);

L = PressureLaplacian_withShape(Mesh,Square);

q = L\ones(Nx*Ny,1);

imagesc(reshape(q,Nx,Ny))

condest(L)




end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end