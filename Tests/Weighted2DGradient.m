function Weighted2DGradient()

%% Compute 1st derivative of u = x^3 +x^2 -1 +y^3 +y^2 -y on non-uniform mesh
%% with 4 point stencil i-1 i i+1 i+2, j-1 j j+1 j+2

%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

Nx = 40;
Ny = 40;
W = 4;
H = 4;
% growth rate
wx = 1.5;
wy = 1.5;

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

Dx = Gradient_Nonuniform(Nx,dx_centre,dx);
Dy = Gradient_Nonuniform(Ny,dy_centre,dy);

Ix = speye(Nx,Nx);
Iy = speye(Ny,Ny);

D_xy = kron(Iy,Dx) + kron(Dy,Ix);

udata = zeros(Nx*Ny,1);
Gradu = zeros(Nx*Ny,1);

for i = 1:length(x_centre)
    for j = 1:length(y_centre)

     udata(index(i,j,Nx)) = x_centre(i)^3 +x_centre(i)^2 -x_centre(i) ...
                            +y_centre(j)^3+y_centre(j)^2 -y_centre(j);
     Gradu(index(i,j,Nx)) = 3*x_centre(i)^2 +2*x_centre(i) -1 ...
                            +3*y_centre(j)^2 +2*y_centre(j) -1;
    end
end

derivative = D_xy * udata;



[Xq,Yq] = meshgrid(W*(1:Nx)./Nx,H*(1:Ny)./Ny);
[X,Y] = meshgrid(x_centre,y_centre);

vq = griddata(X,Y,reshape(derivative,Nx,Ny),Xq,Yq);
Exactq = griddata(X,Y,reshape(Gradu,Nx,Ny),Xq,Yq); 

figure(1)
PlotMesh(Mesh)

figure(2) 
imagesc(H*(1:Ny)./Ny,W*(1:Nx)./Nx,vq)
view(-90,90)
colormap(turbo)
set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

figure(3)
imagesc(H*(1:Ny)./Ny,W*(1:Nx)./Nx,Exactq);
view(-90,90)
colormap(turbo)
set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

figure(4)
imagesc(full(D_xy))
set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
set(gcf,'color','w');




end

function k = index(i,j,Nx)
k = i + (j-1)*Nx;
end
