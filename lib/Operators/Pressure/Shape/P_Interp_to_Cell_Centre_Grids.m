function [p_x_Grid,p_y_Grid,Xq,Yq] = ...
P_Interp_to_Cell_Centre_Grids(Mesh)

Nx = Mesh.Nx;
Ny = Mesh.Ny;

p_x_Grid = zeros(Nx,Ny);
p_y_Grid = zeros(Nx,Ny);

for x=1:Nx
p_y_Grid(x,:) = Mesh.P_centre_y;
end

for y = 1:Ny
p_x_Grid(:,y) = Mesh.P_centre_x;
end

dx_uniform = Mesh.W/Nx;
dy_uniform = Mesh.H/Ny;

[Xq,Yq] = meshgrid((1:Nx).*dx_uniform,(1:Ny).*dy_uniform);

end