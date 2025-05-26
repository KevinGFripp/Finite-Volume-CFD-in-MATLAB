function Ru = Interpolate_u_velocity_to_v_grid(Mesh)
%% Linearly interpolate u grid onto Pressure grid v faces on the interior
%% Pressure grid v faces have known values from the boundary conditions
%% number of columns ((Nx+1) * Ny)
%% number of rows (Nx)*(Ny+1)

Nx = Mesh.Nx;
Ny = Mesh.Ny;

Rx = zeros(Nx,Nx+1);
for i = 1:Nx
 dsi = 1/2 *(Mesh.u_centre_x(i) +Mesh.u_centre_x(i+1)) - Mesh.P_centre_x(i);
 di =  Mesh.u_dx_centre(i);
 coefficientsi = [(1/2 + dsi/di);(1/2 -dsi/di)];
 indsi = i:i+1;
 Rx(i,indsi) = coefficientsi;
end
Rx = sparse(Rx);


Ry = zeros(Ny+1,Ny);
for j = 2:Ny
%    dsj = Mesh.v_centre_y(j) - Mesh.yface(j);
%    dj = Mesh.v_dy(j);
   coefficientsj = [(1/2); (1/2)];
   indsj = j-1:j;
   Ry(j,indsj) = coefficientsj;
end
Ry = sparse(Ry);

Ru = kron(Ry,Rx);




end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end