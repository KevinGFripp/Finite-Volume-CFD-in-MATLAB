function TestCellFaceToCellFaceGradient()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

Width = 1.0;
Height = 1.0;
Nx = 4;
Ny = 4;
GrowthRatesx = 1.0;
GrowthRatesy = 1.0;

Mesh = MakeRectilinearMesh(Nx,Ny,Width,Height,GrowthRatesx,GrowthRatesy,[]);

 figure(3)
 PlotMesh(Mesh)
 hold on
 PlotMeshCentres(Mesh)
 hold off


end