function [v_x_Grid,v_y_Grid,Xq,Yq] = ...
v_velocity_Interp_to_Cell_Centre_Grids(Mesh)

Nx = Mesh.Nx;
Ny = Mesh.Ny;

v_x_Grid = zeros(Nx,Ny+1);
v_y_Grid = zeros(Nx,Ny+1);

for x=1:Nx
  v_y_Grid(x,:) = Mesh.v_centre_y;
end

for y = 1:Ny+1
  v_x_Grid(:,y) = Mesh.v_centre_x;
end

dx_uniform = Mesh.W/Nx;
dy_uniform = Mesh.H/Ny;

[Xq,Yq] = meshgrid((1:Nx).*dx_uniform,(1:Ny).*dy_uniform);

end