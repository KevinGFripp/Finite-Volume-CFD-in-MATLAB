function Rv = Interpolate_v_velocity_to_u_grid(Mesh)
%% Linearly interpolate v grid onto Pressure grid u faces on the interior
%% Pressure grid u faces have known values from the boundary conditions
%% number of columns (Nx * Ny+1)
%% number of rows (Nx+1)*Ny

Nx = Mesh.Nx;
Ny = Mesh.Ny;

Rx = zeros(Nx+1,Nx);
for i = 2:Nx
%   dsi = Mesh.u_centre_x(i) - Mesh.xface(i);
%   di = 1/2 * (Mesh.dx(i-1) + Mesh.dx(i));
%   coefficientsi = [(1/2 + dsi/di); (1/2 - dsi/di)];
  coefficientsi = [1/2;1/2];
  indsi = i-1:i;
  Rx(i,indsi) = coefficientsi;
end
Rx = sparse(Rx);


Ry = zeros(Ny,Ny+1);
for j = 1:Ny
     dsj = 1/2 *(Mesh.v_centre_y(j)+Mesh.v_centre_y(j+1)) ...
           - Mesh.P_centre_y(j);
     dj =  Mesh.v_dy_centre(j);
     coefficientsj = [(1/2 + dsj/dj);(1/2 -dsj/dj)];
     indsj = j:j+1;
     Ry(j,indsj) = coefficientsj;
end
Ry = sparse(Ry);

Rv = kron(Ry,Rx);




end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end