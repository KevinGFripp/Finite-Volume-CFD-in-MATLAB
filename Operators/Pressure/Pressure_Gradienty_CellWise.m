function Dy = Pressure_Gradienty_CellWise(Mesh)

%% 1D Gradient centred-difference, Shifted by one, 
%%  computed on non-uniform cartesian grid
%% Cell face to Cell Face, derivative at cell centre
%% No Boundary Conditions
%% No need to compute coefficients
% dy_centre = distance between cell centroids -> u grid cell size

D = zeros(Mesh.Ny+1,Mesh.Ny);

%left-boundary
% Dy =0

%% Interior
for n = 1:Mesh.Ny-1
  inds = (n):(n+1);
  coeff =[-1;1]./Mesh.P_dy_centre(n);
  D(n+1,inds) = coeff;
end
%%
%right boundary Dy = 0

D = sparse(D);

Ix = speye(Mesh.Nx,Mesh.Nx);

Dy = kron(D,Ix);


end
