function Dx = Gradient_CellWise_Nonuniform(Nx,dx_centre)

%% 1D Gradient centred-difference computed on non-uniform cartesian grid
%% Cell face to Cell Face, derivative at cell centre
%% No Boundary Conditions
%% No need to compute coefficients
% dx_centre = distance between cell centroids

Dx = zeros(Nx-1,Nx-1);

%% Interior
for n = 1:Nx-1

  inds = (n):(n+1);
  coeff =[-1;1]./dx_centre(n);
  Dx(n,inds) = coeff;
end
%%

Dx = sparse(Dx);

end
