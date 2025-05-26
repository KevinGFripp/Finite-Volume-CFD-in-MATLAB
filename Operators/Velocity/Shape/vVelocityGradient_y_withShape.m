function Dv = vVelocityGradient_y_withShape(Mesh,Shape)
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

 Dvy = Gradient_Nonuniform(Mesh.Ny+1, ...
                           Mesh.v_dy_centre, ...
                           Mesh.v_dy);


Ix = speye(Mesh.Nx,Mesh.Nx);
Dv = kron(Dvy,Ix);

%% Shape 
igrid = Shape.v_grid_i;
jgrid = Shape.v_grid_j;

%% Interior set to zero

for i = 1 : length(igrid)
    for j= 2 : (length(jgrid) - 1)
    ind = index(igrid(i),jgrid(j),Mesh.Nx);
    Dv(ind,:) = 0;
    end
end

%% Top Boundary
m = jgrid(end);

for i=1:length(igrid)

   ind = index(igrid(i),m,Mesh.Nx);
   stencil = [(Mesh.v_dy_centre(m-1) + Mesh.v_dy_centre(m-2)) ...
              Mesh.v_dy_centre(m-1)];
   Cn = SolveForCoefficient_RightEdge(stencil,Mesh.v_dy(m));
   inds = index(igrid(i),(m-2):(m),Mesh.Nx);
   coeff =[Cn(1) Cn(2) -sum(Cn)];
   Dv(ind,:) = 0;
   Dv(ind,inds) = coeff;

end


%% Bottom Boundary
m = jgrid(1);

for i=1:length(igrid)
  ind = index(igrid(i),m,Mesh.Nx);
  stencil = [Mesh.v_dy_centre(m) ...
            (Mesh.v_dy_centre(m)+Mesh.v_dy_centre(m+1))];
  Cn = SolveForCoefficient_LeftEdge(stencil,Mesh.v_dy(m));
  inds = index(igrid(i),(m):(m+2),Mesh.Nx);
  coeff =[-sum(Cn) Cn(1) Cn(2)];
  Dv(ind,:) = 0;
  Dv(ind,inds) = coeff;

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
