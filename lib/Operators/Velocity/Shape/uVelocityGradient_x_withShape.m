function Du = uVelocityGradient_x_withShape(Mesh,Shape)
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

%% 
Dux = Gradient_Nonuniform(Mesh.Nx+1, ...
                                  Mesh.u_dx_centre, ...
                                  Mesh.u_dx);


Iy = speye(Mesh.Ny,Mesh.Ny);
Du = kron(Iy,Dux);

%% Shape 
igrid = Shape.u_grid_i;
jgrid = Shape.u_grid_j;

%% Interior set to zero
for i = 2 : (length(igrid) - 1)
    for j= 1 : length(jgrid)
    ind = index(igrid(i),jgrid(j),Mesh.Nx+1);
    Du(ind,:) = 0;
    end
end

%% Left Boundary
m = igrid(1);

for j=1:length(jgrid)
   ind = index(m,jgrid(j),Mesh.Nx+1);
   stencil = [(Mesh.u_dx_centre(m-1) + Mesh.u_dx_centre(m-2)) ...
              Mesh.u_dx_centre(m-1)];
   Cn = SolveForCoefficient_RightEdge(stencil,Mesh.u_dx(m));
   inds = index((m-2):(m),jgrid(j),Mesh.Nx+1);
   coeff =[Cn(1) Cn(2) -sum(Cn)];
   Du(ind,:) = 0;
   Du(ind,inds) = coeff;

end


%% Right Boundary
m = igrid(end);

for j=1:length(jgrid)
  ind = index(m,jgrid(j),Mesh.Nx+1);
  stencil = [Mesh.u_dx_centre(m) ...
            (Mesh.u_dx_centre(m)+Mesh.u_dx_centre(m+1))];
  Cn = SolveForCoefficient_LeftEdge(stencil,Mesh.u_dx(m));
  inds = index((m):(m+2),jgrid(j),Mesh.Nx+1);
  coeff =[-sum(Cn) Cn(1) Cn(2)];
  Du(ind,:) = 0;
  Du(ind,inds) = coeff;

end


end

function C = SolveForCoefficient_LeftEdge(stencil,dx)
% %% 3 point stencil, 2nd order accurate
% % ah0 +bh1  = 2*dx
% % +ah0^2 +bh1^2 = 0
CoefficientMatrix = zeros(2,2);
 b = zeros(2,1);
 b(1,1) = 2*dx;

 CoefficientMatrix(1,:) =[stencil(1) stencil(2)];
 CoefficientMatrix(2,:) =[stencil(1) stencil(2)].^2;
 
  C = CoefficientMatrix\b;
  C = C./2./dx;
end

function C = SolveForCoefficient_RightEdge(stencil,dx)
% %% 3 point stencil, 2nd order accurate
% % -ah0 -bh1  = 2*dx
% % +ah0^2 +bh1^2 = 0
CoefficientMatrix = zeros(2,2);
 b = zeros(2,1);
 b(1,1) = 2*dx;

 CoefficientMatrix(1,:) =[-stencil(1) -stencil(2)];
 CoefficientMatrix(2,:) =[stencil(1) stencil(2)].^2;
 
  C = CoefficientMatrix\b;
  C = C./2./dx;
end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end