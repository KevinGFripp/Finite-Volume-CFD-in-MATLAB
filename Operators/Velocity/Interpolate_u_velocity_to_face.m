function Ru = Interpolate_u_velocity_to_face(Mesh)
Nx = Mesh.Nx;
Ny = Mesh.Ny;

Rx = zeros(Nx+1,Nx+1);
Rx(1,1)=1;
Rx(Nx+1,Nx+1)=1;
Iy = speye(Ny,Ny);

% interior points need to be interpolated onto the Pressure grid faces
for n =2:Nx

  Offset = (Mesh.u_centre_x(n) - Mesh.xface(n));

 if(Offset >= 0.0)
  inds = n-1:n;
  separation = 1/2 * (Mesh.u_dx(n-1) + Mesh.u_dx(n));
  coeff1 = (Offset/separation);
  coeff2 = (separation - Offset)/separation;
  coefficients = [coeff1;coeff2];
 else
 inds = n:n+1;
  separation = 1/2 * (Mesh.u_dx(n) + Mesh.u_dx(n+1));
  coeff1 = (abs(Offset)/separation);
  coeff2 = (separation - abs(Offset))/separation;
  coefficients = [coeff2;coeff1];
 end

Rx(n,inds) = coefficients;
end

Rx = sparse(Rx);

Ru = kron(Iy,Rx);



end