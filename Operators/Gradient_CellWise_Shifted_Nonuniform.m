function Dx = Gradient_CellWise_Shifted_Nonuniform(Nx,dx_centre)

%% 1D Gradient centred-difference, Shifted by one, 
%%  computed on non-uniform cartesian grid
%% Cell face to Cell Face, derivative at cell centre
%% No Boundary Conditions
%% No need to compute coefficients
% dx_centre = distance between cell centroids

Dx = zeros(Nx-1,Nx-1);

%% Interior
for n = 2:Nx

  inds = (n-1):(n);
  coeff =[-1;1]./dx_centre(n);
  Dx(n,inds) = coeff;
end
%%

Dx = sparse(Dx);

end
