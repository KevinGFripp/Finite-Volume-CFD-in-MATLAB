function Dx = Gradient_Nonuniform_Neumann(Nx,dx_centre,dx)

%% 1D Gradient computed on non-uniform cartesian grid
%% Dirichlet Boundary Conditions Everywhere
% dx_centre = distance between cell centroids in x or y direction
% dx = distance between cell faces in x or y direction

Dx = zeros(Nx,Nx);


%% Left Boundary
  m = 1;
  stencil = [dx(m) dx_centre(m)];
  Cn = SolveForCoefficient(stencil,dx(m));
  inds = (m):(m+1);
  coeff =[-(sum(Cn)-Cn(1));Cn(2)];
  Dx(m,inds) = coeff;
%%


%% Interior
for n = 2:Nx-1
  stencil = [dx_centre(n-1) dx_centre(n)];
  Cn = SolveForCoefficient(stencil,dx(n));
  inds = (n-1):(n+1);
  coeff =[Cn(1);-sum(Cn);Cn(2)];
  Dx(n,inds) = coeff;
end
%%


%% Right Boundary
  m = Nx;
  stencil = [dx_centre(m-1) dx(m)];
  Cn = SolveForCoefficient(stencil,dx(m));
  inds = (m-1):(m);
  coeff =[Cn(1);-(sum(Cn)-Cn(2))];
  Dx(m,inds) = coeff;
%%


Dx = sparse(Dx);

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
