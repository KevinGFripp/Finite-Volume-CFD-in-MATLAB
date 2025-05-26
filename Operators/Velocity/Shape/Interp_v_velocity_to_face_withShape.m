function Rv = Interp_v_velocity_to_face_withShape(Mesh,Shape)

Nx = Mesh.Nx;
Ny = Mesh.Ny;

Ry = zeros(Ny+1,Ny+1);
Ry(1,1)=1;
Ry(Ny+1,Ny+1)=1;
Ix = speye(Nx,Nx);


% interior points need to be interpolated onto the Pressure grid faces
for n =2:Ny

 Offset = (Mesh.v_centre_y(n) - Mesh.yface(n));

 if(Offset >= 0.0)
  inds = n-1:n;
  separation = 1/2 * (Mesh.v_dy(n-1) + Mesh.v_dy(n));
  coeff1 = (Offset/separation);
  coeff2 = (separation - Offset)/separation;
  coefficients = [coeff1;coeff2];
 else
 inds = n:n+1;
  separation = 1/2 * (Mesh.v_dy(n) + Mesh.v_dy(n+1));
  coeff1 = (abs(Offset)/separation);
  coeff2 = (separation - abs(Offset))/separation;
  coefficients = [coeff2;coeff1];
 end
 
Ry(n,inds) = coefficients;
end

Ry = sparse(Ry);
Rv = kron(Ry,Ix);



end