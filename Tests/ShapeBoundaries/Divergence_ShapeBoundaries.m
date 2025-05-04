function Divergence_ShapeBoundaries()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

% Grid
Nx = 16;
Ny = 16;
Width = 1;
Height = 1;
% growth rate
wx = 0.00008;
wy = 0.00008;

Square = DefineBody(0.25,0.25,Width/Nx,Height/Ny,Nx,Ny,Width/Nx,Height/Ny);
Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height,[wx wx wx],[wy wy wy],Square);

% PlotMeshandShape(Mesh,Square)
Dx = Divergence_x_withShape(Mesh,Square);
Dy = Divergence_y_withShape(Mesh,Square);

Datax = linspace(0,1,(Nx+1)*Ny).';
Datay = linspace(0,1,Nx*(Ny+1)).';

subplot(1,3,1)
imagesc(reshape(Dx*Datax,Nx,Ny))
view(-90,90)

set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

subplot(1,3,2)
imagesc(reshape(Dy*Datay,Nx,Ny))
view(-90,90)

set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

set(gcf,'Position',[800 500 700 300]);
colormap(jet)

subplot(1,3,3)
imagesc(reshape(Datax,Nx+1,Ny))
view(-90,90)

end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end