function u_Gradient_y_ShapeBoundaries()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

% Grid
Nx = 10;
Ny = 10;
Width = 1;
Height = 1;
% growth rate
wx = 0.0005;
wy = 0.0005;

Square = DefineBody(0.25,0.25,Width/Nx,Height/Ny,Nx,Ny,Width/Nx,Height/Ny)
Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height,[wx wx wx],[wy wy wy],Square);

Boundaries = CreateBoundaries('Wall',0,0,'Wall',0,0,'Wall',0,0,'Wall',0,0);

Duy = uVelocityGradient_y_withShape(Mesh,Boundaries,Square);

Datax = zeros((Nx+1)*Ny,1);

for i=1:Nx+1
    for j =1:Ny
    Datax(index(i,j,Nx+1)) = (1);  
    end
end

figure(2)
imagesc(reshape(Duy*Datax,Nx+1,Ny))
view(-90,90)
set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

colormap(jet)

end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end