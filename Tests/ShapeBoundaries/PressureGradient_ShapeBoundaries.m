function PressureGradient_ShapeBoundaries()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

% Grid
Nx = 40;
Ny = 40;
Width = 1;
Height = 1;
% growth rate
wx = 0.8;
wy = 0.8;

Square = DefineBody(0.25,0.25,Width/Nx,Height/Ny,Nx,Ny,Width/Nx,Height/Ny);
Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height,[wx wx wx],[wy wy wy],Square);

% PlotMeshandShape(Mesh,Square)
Dx = Pressure_Gradientx_CellWise_withShape(Mesh,Square);
Dy = Pressure_Gradienty_CellWise_withShape(Mesh,Square);

Datax = linspace(0,1,Nx*Ny).';
Datay = linspace(0,1,Nx*Ny).';

subplot(1,2,1)
imagesc(reshape(Dx*Datax,Nx+1,Ny))
view(-90,90)

set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

subplot(1,2,2)
imagesc(reshape(Dy*Datay,Nx,Ny+1))
view(-90,90)

set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

set(gcf,'Position',[800 500 700 300]);
colormap(jet)

end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end