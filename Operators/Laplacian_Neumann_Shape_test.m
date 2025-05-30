function Lx = Laplacian_Neumann_Shape_test(Nx,dx_centre,dx,Shape_Grid_Indexes)

%% 1D Laplacian computed on non-uniform cartesian grid
%% Neumann Boundary Conditions everywhere
%% Deals with Boundary of Shape -> Shape_Grid_Indexes
%% Shape_Grid_Indexes(1) = start of boundary from the left
%% Shape_Grid_Indexes(end) = end of boundary from the left
%% interior Shape_Grid_Indexes(2:end-1) are the identity matrix
% dx_centre = distance between cell centroids in x or y direction
% dx = distance between cell faces in x or y direction

Lx = zeros(Nx,Nx);

% %% Interior
% for n=2:Nx-2  
%   stencil = [dx_centre(n-1) dx_centre(n) (dx_centre(n) + dx_centre(n+1))];
%   Cn = SolveForCoefficient(stencil,dx(n));
%   coeff = [Cn(1);-sum(Cn);Cn(2);Cn(3)];
%   inds = (n-1):(n+2);
% Lx(n,inds) = coeff;
% end
% %%
% 
% %% Left Boundary
% m = 1;
% stencil = [dx(m) dx_centre(m) (dx_centre(m) +dx_centre(m+1))];
% Cm = SolveForCoefficient(stencil,dx(m));
% coeff = [-(sum(Cm)-Cm(1));Cm(2);Cm(3)];
% inds = (m):(m+2);
% Lx(m,inds) = coeff;
% %%
% 
% 
% %% One-from Right Boundary
% m = Nx-1;
% stencil = [(dx_centre(m-2) +dx_centre(m-1)) ...
%            (dx_centre(m-1)) ...
%             dx_centre(m)];
% Cm = SolveForCoefficient_OneFromRightEdge(stencil,dx(m));
% coeff = [Cm(1);Cm(2);-sum(Cm);Cm(3)];
% inds = (m-2):(m+1);
% Lx(m,inds) = coeff;
% %%
% 
% 
% %% Right Boundary
% m = Nx;
% stencil = [(dx_centre(m-2) +dx_centre(m-1)) dx_centre(m-1) dx(m)];
% Cm = SolveForCoefficient_OneFromRightEdge(stencil,dx(m));
% coeff = [Cm(1);Cm(2);-(sum(Cm)-Cm(3))];
% inds = (m-2):(m);
% Lx(m,inds) = coeff;
% %%

%% Shape Boundaries

% Interior
Shape_index_interior = Shape_Grid_Indexes(1:end);

for s = 1:length(Shape_index_interior)

    ind_diagonal = Shape_index_interior(s);
    Lx(ind_diagonal,:) = 0;
    Lx(ind_diagonal,ind_diagonal) = 1;

end

% % left wall boundary
 m = Shape_Grid_Indexes(1)-1;
 stencil = [(dx_centre(m-2) +dx_centre(m-1)) dx_centre(m-1) dx(m)];
 Cm = SolveForCoefficient_OneFromRightEdge(stencil,dx(m));
 coeff = [1;1;1];
 inds = (m-2):(m);
 Lx(m,:) = 0;
 Lx(m,inds) = coeff;
 
 
 % one from left boundary
 m = Shape_Grid_Indexes(1)-2;
 stencil = [(dx_centre(m-2) +dx_centre(m-1)) ...
           (dx_centre(m-1)) ...
            dx_centre(m)];
 Cm = SolveForCoefficient_OneFromRightEdge(stencil,dx(m));
 coeff = [1;1;1;1];
 inds = (m-2):(m+1);
 Lx(m,:) = 0;
 Lx(m,inds) = coeff;


% % right boundary
m = Shape_Grid_Indexes(end)+1;
stencil = [dx(m) dx_centre(m) (dx_centre(m) +dx_centre(m+1))];
Cm = SolveForCoefficient(stencil,dx(m));
coeff = [1;1;1];
inds = (m):(m+2);
Lx(m,:) = 0;
Lx(m,inds) = coeff;


Lx = sparse(Lx);

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

function C = SolveForCoefficient_LeftEdge(stencil,dx)

% %% 3 point stencil, 1st order accurate
% % bh1 +ch2 = 0
% % bh1^2 +ch2^2 = 2dx1^2

CoefficientMatrix = zeros(2,2);
 b = zeros(2,1);
 b(2,1) = 2 * dx.^2;

 CoefficientMatrix(1,:) =[stencil(1) stencil(2)];
 CoefficientMatrix(2,:) =[stencil(1) stencil(2)].^2;

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

function C = SolveForCoefficient_RightEdge(stencil,dx)

% %% 3 point stencil, 1rd order accurate
% % -ah0 -bh1 = 0
% % +ah0^2 +bh1^2 = 2h0^2

CoefficientMatrix = zeros(2,2);
 b = zeros(2,1);
 b(2,1) = 2 * dx.^2;

 CoefficientMatrix(1,:) =[-stencil(1) -stencil(2)];
 CoefficientMatrix(2,:) =[stencil(1) stencil(2)].^2;

  C = CoefficientMatrix\b;
  C = C./dx./dx;
end
