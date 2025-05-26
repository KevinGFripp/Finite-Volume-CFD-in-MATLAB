function [u_x_Grid,u_y_Grid,Xq,Yq] = ...
u_velocity_Interp_to_Cell_Centre_Grids(Mesh)

Nx = Mesh.Nx;
Ny = Mesh.Ny;

u_x_Grid = zeros(Nx+1,Ny);
u_y_Grid = zeros(Nx+1,Ny);

for x=1:Nx+1
u_y_Grid(x,:) = Mesh.u_centre_y;
end

for y = 1:Ny
u_x_Grid(:,y) = Mesh.u_centre_x;
end

dx_uniform = Mesh.W/Nx;
dy_uniform = Mesh.H/Ny;

[Xq,Yq] = meshgrid((1:Nx).*dx_uniform,(1:Ny).*dy_uniform);

end