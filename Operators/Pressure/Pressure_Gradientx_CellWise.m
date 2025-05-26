function Dx = Pressure_Gradientx_CellWise(Mesh)

%% 1D Gradient centred-difference, Shifted by one, 
%%  computed on non-uniform cartesian grid
%% Cell face to Cell Face, derivative at cell centre
%% No Boundary Conditions
%% No need to compute coefficients
% dx_centre = distance between cell centroids -> u grid cell size

D = zeros(Mesh.Nx+1,Mesh.Nx);

%left-boundary
% Dx =0

%% Interior
for n = 1:Mesh.Nx-1
  inds = (n):(n+1);
  coeff =[-1;1]./Mesh.P_dx_centre(n);
  D(n+1,inds) = coeff;
end
%%
%right boundary Dx = 0

D = sparse(D);

Iy = speye(Mesh.Ny,Mesh.Ny);

Dx = kron(Iy,D);


end
