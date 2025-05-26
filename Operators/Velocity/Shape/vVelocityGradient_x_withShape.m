function [Dv,bv_bc] = vVelocityGradient_x_withShape(Mesh,Boundaries,Shape)
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%


%% velocity boundary conditions on left and right
[Dvx,bx_bc] = Gradient_Nonuniform_Dirichlet(Mesh.Nx, ...
                                     Mesh.v_dx_centre, ...
                                     Mesh.v_dx,...
                                     Boundaries.Left_v, ...
                                     Boundaries.Right_v);

if(strcmp(Boundaries.RightWallType,'Free'))
Dvx_Nm = Gradient_Nonuniform_Neumann(Mesh.Nx, ...
                                     Mesh.v_dx_centre, ...
                                     Mesh.v_dx);
 Dvx(Mesh.Nx,:) = Dvx_Nm(Mesh.Nx,:);
end

Iy = speye(Mesh.Ny+1,Mesh.Ny+1);

Dv = kron(Iy,Dvx);

bv_bc = zeros(Mesh.Nx*(Mesh.Ny+1),1);

LeftWall = index(1:Mesh.Nx,1,Mesh.Nx);
RightWall = index(1:Mesh.Nx,Mesh.Ny+1,Mesh.Nx);

for n= 1:length(LeftWall)
   bv_bc(LeftWall(n)) = bx_bc(1);
   bv_bc(RightWall(n)) = bx_bc(Mesh.Nx);
end

%% Shape 
igrid = Shape.v_grid_i;
jgrid = Shape.v_grid_j;

%% interior set to 0
for i = 1:length(igrid)
    for j=1:length(jgrid)    
        ind = index(igrid(i),jgrid(j),Mesh.Nx);
        Dv(ind,:) = 0;
    end
end


if(strcmp(Boundaries.RightWallType,'Wall') ...
   || strcmp(Boundaries.RightWallType,'Inlet'))
%% Right wall 0 Dirichlet boundary
m = igrid(end)+1;

for j =1:length(jgrid)
  ind = index(m,jgrid(j),Mesh.Nx);  
  stencil = [Mesh.v_dx(m) Mesh.v_dx_centre(m)];
  Cn = SolveForCoefficient(stencil,Mesh.v_dx(m));
  inds = index((m):(m+1),jgrid(j),Mesh.Nx);
  coeff =[-(sum(Cn)+Cn(1)) Cn(2)];
  Dv(ind,:) = 0;
  Dv(ind,inds) = coeff;
end

end

if(strcmp(Boundaries.LeftWallType,'Wall') ...
   || strcmp(Boundaries.LeftWallType,'Inlet'))
%% Left wall 0 Dirichlet boundary
m = igrid(1)-1;

for j =1:length(jgrid)
  ind = index(m,jgrid(j),Mesh.Nx);  
  stencil = [Mesh.v_dx_centre(m-1) Mesh.v_dx(m)];
  Cn = SolveForCoefficient(stencil,Mesh.v_dx(m));
  inds = index((m-1):(m),jgrid(j),Mesh.Nx);
  coeff =[Cn(1) -(sum(Cn)+Cn(2))];
  Dv(ind,:) = 0;
  Dv(ind,inds) = coeff;
end

end



end

function C = SolveForCoefficient(stencil,dx)

% %% 3 point stencil, 2nd order accurate
% % -ah0 +bh1  = 2* dx
% % +ah0^2 +bh1^2 = 0

CoefficientMatrix = zeros(2,2);
 b = zeros(2,1);
 b(1,1) = 2*dx;

 CoefficientMatrix(1,:) =[-stencil(1) stencil(2)];
 CoefficientMatrix(2,:) =[stencil(1) stencil(2)].^2;
 
  C = CoefficientMatrix\b;
  C = C./2./dx;
end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end