function PoissonSolution_Boundaries()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

% Grid
Nx = 80;
Ny = 80;
W = 1;
H = 1;
% growth rate
wx = 1;
wy = 1;
Mesh = MakeRectilinearMesh(Nx,Ny,W,H,wx,wy,[]);

% Mesh -> non-uniform position of cell faces
xn = Mesh.xface;
yn = Mesh.yface;

% cell centroids
x_centre = (xn(2:Nx+1)+xn(1:Nx))./2;
y_centre = (yn(2:Ny+1)+yn(1:Ny))./2;

% cell distance between cell centres
dx_centre = x_centre(2:Nx) - x_centre(1:Nx-1);
dy_centre = y_centre(2:Ny) - y_centre(1:Ny-1);

% grid cell size
dx = (xn(2:Nx+1) - xn(1:Nx));
dy = (yn(2:Ny+1) - yn(1:Ny));

disp(strcat('Minimum Cell Size = ',num2str(min(min(dx,[],'all'),min(dy,[],'all')))));
disp(strcat('Maximum Cell Size = ',num2str(max(max(dx,[],'all'),max(dy,[],'all')))));


% Dirichlet boundary conditions on all sides
Lx = Laplacian_Nonuniform_Dirichlet(Nx,dx_centre,dx,0,0);
[Ly,by_bc] = Laplacian_Nonuniform_Dirichlet(Ny,dy_centre,dy,0,-1.0);



Ix = speye(Nx,Nx);
Iy = speye(Ny,Ny);
L_xy = kron(Iy,Lx) + kron(Ly,Ix);
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
PlotMesh(Mesh)



end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end

