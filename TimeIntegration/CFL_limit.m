function maxdt = CFL_limit(LengthScale,Mesh,Parameters)
% Nu * dt / dx**2 <= 1 CFL condition

maxdt = ((LengthScale).^2) * ...
      min( (min(Mesh.u_dx)^2)/Parameters.Nu , ...
         (min(Mesh.u_dy)^2)/Parameters.Nu ); 

end