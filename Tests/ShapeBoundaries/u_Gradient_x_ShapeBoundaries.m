function u_Gradient_x_ShapeBoundaries()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

% Grid
Nx = 32;
Ny = 32;
Width = 1;
Height = 1;
% growth rate
wx = 0.0005;
wy = 0.0005;

Square = DefineBody(0.25,0.25,Width/Nx,Height/Ny,Nx,Ny,Width/Nx,Height/Ny);
Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height,[wx wx wx],[wy wy wy],Square);

% PlotMeshandShape(Mesh,Square)
% hold on
% PlotMeshCentres(Mesh)

Dux = uVelocityGradient_x_withShape(Mesh,Square);

Datax = linspace(0,1,Nx*(Ny+1)).';


figure(2)
imagesc(reshape(Dux*Datax,Nx+1,Ny))
view(-90,90)
set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

colormap(jet)

end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end