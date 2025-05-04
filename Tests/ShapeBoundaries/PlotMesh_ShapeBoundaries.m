function PlotMesh_ShapeBoundaries()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

% Grid
Nx = 6;
Ny = 6;
Width = 1;
Height = 1;
% growth rate
wx = 0.8;
wy = 0.8;

Square = DefineBody(0.25,0.25,Width/Nx,Height/Ny,Nx,Ny,Width/Nx,Height/Ny);
Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height,[wx wx wx],[wy wy wy],Square);
PlotMeshandShape(Mesh,Square)
hold on
PlotMeshCentres(Mesh)
hold off

end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end