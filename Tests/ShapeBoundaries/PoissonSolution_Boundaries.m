function PoissonSolution_Boundaries()
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
wx = 0.0001;
wy = 0.0001;

Square = DefineBody(0.1,0.1,Width/Nx,Height/Ny,Nx,Ny,Width/Nx,Height/Ny);

Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height,[wx wx wx],[wy wy wy],Square);

Lx_s = Laplacian_Neumann_Shape(Nx,Mesh.P_dx_centre,Mesh.dx,Square.p_grid_i);
Ly_s = Laplacian_Neumann_Shape(Ny,Mesh.P_dy_centre,Mesh.dy,Square.p_grid_j);

% Dirichlet boundary conditions on all sides
Lx = Laplacian_Nonuniform_Dirichlet(Nx,Mesh.P_dx_centre,Mesh.dx,0,0);
[Ly,by_bc] = Laplacian_Nonuniform_Dirichlet(Ny,Mesh.P_dy_centre,Mesh.dy,0,-1.0);

Lx = Lx_s;
Ly= Ly_s;

Ix = speye(Nx,Nx);
Iy = speye(Ny,Ny);
L_xy = kron(Iy,Lx) + kron(Ly,Ix);
L_xy(1,:) = 0;
L_xy(1,1) = 1;
disp(strcat('Condition number = ',num2str(condest(L_xy))))


b = zeros(Nx*Ny,1);
for x = 1:Nx
b(index(x,Ny,Nx)) = by_bc(Ny);
end


t=L_xy\b;

ax1 = figure(1);
imagesc(reshape(t,Nx,Ny));
set(gca,'FontSize',16,'FontName','Times','fontweight','normal');
set(gcf,'color','w');
view(-90,90)
colormap(ax1,turbo);

figure(2)
imagesc(full(L_xy))




end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end

