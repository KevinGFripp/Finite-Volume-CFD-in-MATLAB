function v_Gradient_x_ShapeBoundaries()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

% Grid
Nx = 24;
Ny = 24;
Width = 1;
Height = 1;
% growth rate
wx = 0.75;
wy = 0.75;

Square = DefineBody(0.25,0.25,Width/Nx,Height/Ny,Nx,Ny,Width/Nx,Height/Ny);
Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height,[wx wx wx],[wy wy wy],Square);

Boundaries = CreateBoundaries('Wall',0,0,'Wall',0,0,'Wall',0,0,'Wall',0,0);

[Dvx,~] = vVelocityGradient_x_withShape(Mesh,Boundaries,Square);

Datax = zeros(Nx*(Ny+1),1);

for i=1:Nx
    for j =1:Ny+1
    Datax(index(i,j,Nx)) = 1;  
    end
end


imagesc(reshape(Dvx*Datax,Nx,Ny+1));
view(-90,90)

set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

colormap(jet)

end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end