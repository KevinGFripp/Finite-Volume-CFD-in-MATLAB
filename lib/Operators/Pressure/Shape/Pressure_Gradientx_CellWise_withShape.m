function Dx = Pressure_Gradientx_CellWise_withShape(Mesh,Shape)

%% 1D Gradient centred-difference, Shifted by one, 
%%  computed on non-uniform cartesian grid
%% Cell face to Cell Face, derivative at cell centre
%% Neumann boudary conditions on shape
% Interior is 0

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

%% Shape interior
% Zero rows
igrid = Shape.u_grid_i;
jgrid = Shape.p_grid_j;
 
for i=1:length(igrid)
    for j=1:length(jgrid)
    ind = index(igrid(i),jgrid(j),Mesh.Nx+1);
    Dx(ind,:) = 0;
    end
end


end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end
