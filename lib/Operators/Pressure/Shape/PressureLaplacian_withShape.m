function L = PressureLaplacian_withShape(Mesh,Shape)
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

%% Pressure Laplacian with Neumann Boundary Conditions on the walls,
%% and on the boundary of the shape

Lx = Laplacian_Nonuniform_Neumann(Mesh.Nx,Mesh.P_dx_centre,Mesh.dx);
Ly = Laplacian_Nonuniform_Neumann(Mesh.Ny,Mesh.P_dy_centre,Mesh.dy);

Iy = speye(Mesh.Ny,Mesh.Ny);
Ix = speye(Mesh.Nx,Mesh.Nx);

Maskxy = speye(Mesh.Nx*Mesh.Ny,Mesh.Nx*Mesh.Ny);


gradientx = kron(Iy,Lx);
gradienty = kron(Ly,Ix);

L = gradientx + gradienty;

%% Shape interior
% diagonal of ones, not solved for
igrid = Shape.p_grid_i;
jgrid = Shape.p_grid_j;

for x = 1:length(igrid)
    for y = 1:length(jgrid)
    
        ind = index(igrid(x),jgrid(y),Mesh.Nx);
        Maskxy(ind,ind) =0;
       
    end
end
%%
L = Maskxy*L;

for x = 1:length(igrid)
    for y = 1:length(jgrid)
    
        ind = index(igrid(x),jgrid(y),Mesh.Nx);
        L(ind,ind) =1;
       
    end
end



%% Shape left edge
% derivative in y is kept
% derivative in x has boundary

% Left Boundary
xval = igrid(1)-1;
for y = 1:length(jgrid)

  ind = index(xval,jgrid(y),Mesh.Nx);
  L(ind,:) = gradienty(ind,:);

 m = xval;
 stencil = [(Mesh.P_dx_centre(m-2) +Mesh.P_dx_centre(m-1)) ...
             Mesh.P_dx_centre(m-1) ...
             Mesh.dx(m)];
 Cm = SolveForCoefficient_OneFromRightEdge(stencil,Mesh.dx(m));
 coeff = [Cm(1) Cm(2) -(sum(Cm)-Cm(3))];
 %inds = (m-2):(m);
 inds = index((xval-2):xval,jgrid(y),Mesh.Nx);
 L(ind,inds) = L(ind,inds) + coeff;

end

% One from left boundary
xval = igrid(1)-2;
for y = 1:length(jgrid)

  ind = index(xval,jgrid(y),Mesh.Nx);
  L(ind,:) = gradienty(ind,:);

 m = xval;
 stencil = [(Mesh.P_dx_centre(m-2) +Mesh.P_dx_centre(m-1)) ...
             (Mesh.P_dx_centre(m-1)) ...
              Mesh.P_dx_centre(m)];
 Cm = SolveForCoefficient_OneFromRightEdge(stencil,Mesh.dx(m));
 coeff = [Cm(1) Cm(2) -sum(Cm) Cm(3)];
 %inds = (m-2):(m+1);
 inds = index((xval-2):(xval+1),jgrid(y),Mesh.Nx);
 L(ind,inds) = L(ind,inds) + coeff;

end
%%





%% Shape right edge
% derivative in y is kept
% derivative in x has boundary

% Right boundary
xval = igrid(end)+1;
for y = 1:length(jgrid)

  ind = index(xval,jgrid(y),Mesh.Nx);
  L(ind,:) = gradienty(ind,:);

m = xval;
stencil = [Mesh.dx(m) ...
           Mesh.P_dx_centre(m) ...
          (Mesh.P_dx_centre(m) +Mesh.P_dx_centre(m+1))];
Cm = SolveForCoefficient(stencil,Mesh.dx(m));
coeff = [-(sum(Cm)-Cm(1)) Cm(2) Cm(3)];
% inds = (m):(m+2);
inds = index(xval:(xval+2),jgrid(y),Mesh.Nx);
L(ind,inds) = L(ind,inds) + coeff;

end
%%





%% Shape bottom edge
% derivative in x is kept
% derivative in y has boundary

% Bottom boundary
yval = jgrid(1)-1;
for x = 1:length(igrid)

  ind = index(igrid(x),yval,Mesh.Nx);
  L(ind,:) = gradientx(ind,:);

 m = yval;
 stencil = [(Mesh.P_dy_centre(m-2) +Mesh.P_dy_centre(m-1)) ...
             Mesh.P_dy_centre(m-1) ...
             Mesh.dy(m)];
 Cm = SolveForCoefficient_OneFromRightEdge(stencil,Mesh.dy(m));
 coeff = [Cm(1) Cm(2) -(sum(Cm)-Cm(3))];
 %inds = (m-2):(m);
 inds = index(igrid(x),(yval-2):yval,Mesh.Nx);
 L(ind,inds) = L(ind,inds) + coeff;

end

% One from bottom boundary
yval = jgrid(1)-2;
for x = 1:length(igrid)

  ind = index(igrid(x),yval,Mesh.Nx);
  L(ind,:) = gradientx(ind,:);

 m=yval;
 stencil = [(Mesh.P_dy_centre(m-2) +Mesh.P_dy_centre(m-1)) ...
             (Mesh.P_dy_centre(m-1)) ...
              Mesh.P_dy_centre(m)];
 Cm = SolveForCoefficient_OneFromRightEdge(stencil,Mesh.dy(m));
 coeff = [Cm(1) Cm(2) -sum(Cm) Cm(3)];
 %inds = (m-2):(m+1);
 inds = index(igrid(x),(yval-2):(yval+1),Mesh.Nx);
 L(ind,inds) = L(ind,inds) + coeff;

end

%%




%% Shape top edge
% derivative in x is kept
% derivative in y has boundary

yval = jgrid(end)+1;
for x = 1:length(igrid)

  ind = index(igrid(x),yval,Mesh.Nx);
  L(ind,:) = gradientx(ind,:);

 m = yval;
stencil = [Mesh.dy(m) ...
           Mesh.P_dy_centre(m) ...
          (Mesh.P_dy_centre(m) +Mesh.P_dy_centre(m+1))];
Cm = SolveForCoefficient(stencil,Mesh.dy(m));
coeff = [-(sum(Cm)-Cm(1)) Cm(2) Cm(3)];
%inds = (m):(m+2);
inds = index(igrid(x),yval:(yval+2),Mesh.Nx);
L(ind,inds) = L(ind,inds) + coeff;

end
%%





L(1,:) = 0;
L(1,1) = 1;
end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end

function C = SolveForCoefficient(stencil,dx)
% %% 4 point stencil, 2nd order accurate
% % -ah0 +bh1 +ch2 = 0
% % -ah0^3 +bh1^3 +ch2^3 = 0
% % +ah0^2 +bh1^2 +ch2^2 = 2h0^2

 CoefficientMatrix = zeros(3,3);
 b = zeros(3,1);
 b(3,1) = 2 * dx.^2;

 CoefficientMatrix(1,:) =[-stencil(1) stencil(2) stencil(3)];
 CoefficientMatrix(2,:) =[-stencil(1) stencil(2) stencil(3)].^3;
 CoefficientMatrix(3,:) =[stencil(1) stencil(2) stencil(3)].^2;
 
  C = CoefficientMatrix\b;
  C = C./dx./dx;
end
function C = SolveForCoefficient_OneFromRightEdge(stencil,dx)
% %% 4 point stencil, 2rd order accurate
% % -ah0 -bh1 +ch2 = 0
% % -ah0^3 -bh1^3 +ch2^3 = 0
% % +ah0^2 +bh1^2 +ch2^2 = 2h0^2

 CoefficientMatrix = zeros(3,3);
 b = zeros(3,1);
 b(3,1) = 2 * dx.^2;

 CoefficientMatrix(1,:) =[-stencil(1) -stencil(2) stencil(3)];
 CoefficientMatrix(2,:) =[-stencil(1) -stencil(2) stencil(3)].^3;
 CoefficientMatrix(3,:) =[stencil(1) stencil(2) stencil(3)].^2;

  C = CoefficientMatrix\b;
  C = C./dx./dx;
end
