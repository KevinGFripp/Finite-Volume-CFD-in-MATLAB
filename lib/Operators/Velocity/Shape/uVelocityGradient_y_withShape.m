function [Du,bu_bc] = uVelocityGradient_y_withShape(Mesh,Boundaries,Shape)

%% velocity boundary conditions on bottom and top
[Duy,by_bc] = Gradient_Nonuniform_Dirichlet(Mesh.Ny, ...
                                     Mesh.u_dy_centre, ...
                                     Mesh.u_dy,...
                                     Boundaries.Bottom_u, ...
                                     Boundaries.Top_u);

Ix = speye(Mesh.Nx+1,Mesh.Nx+1);
Du = kron(Duy,Ix);

bu_bc = zeros((Mesh.Nx+1)*Mesh.Ny,1);

BottomWall = index(1,1:Mesh.Ny,Mesh.Nx+1);
TopWall = index(Mesh.Nx+1,1:Mesh.Ny,Mesh.Nx+1);

for n= 1:length(BottomWall)
   bu_bc(BottomWall(n)) = by_bc(1);
   bu_bc(TopWall(n)) = by_bc(Mesh.Ny);
end

%% Shape 
igrid = Shape.u_grid_i;
jgrid = Shape.u_grid_j;

%% interior set to 0
for i = 1:length(igrid)
    for j=1:length(jgrid)    
        ind = index(igrid(i),jgrid(j),Mesh.Nx+1);
        Du(ind,:) = 0;
    end
end
% 
%% Top wall 0 Dirichlet boundary
m = jgrid(end)+1;

for i =1:length(igrid)
  ind = index(igrid(i),m,Mesh.Nx+1);  
  stencil = [Mesh.u_dy(m) Mesh.u_dy_centre(m)];
  Cn = SolveForCoefficient(stencil,Mesh.u_dy(m));
  inds = index(igrid(i),(m):(m+1),Mesh.Nx+1);
  coeff =[-(sum(Cn)+Cn(1)) Cn(2)];
  Du(ind,:) = 0;
  Du(ind,inds) = coeff;
end

%% Bottom wall 0 Dirichlet boundary
m = jgrid(1)-1;

for i =1:length(igrid)
  ind = index(igrid(i),m,Mesh.Nx+1);  
  stencil = [Mesh.u_dy_centre(m-1) Mesh.u_dy(m)];
  Cn = SolveForCoefficient(stencil,Mesh.u_dy(m));
  inds = index(igrid(i),(m-1):(m),Mesh.Nx+1);
  coeff =[Cn(1) -(sum(Cn)+Cn(2))];
  Du(ind,:) = 0;
  Du(ind,inds) = coeff;
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
