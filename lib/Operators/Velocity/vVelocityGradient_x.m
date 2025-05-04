function [Dv,bv_bc] = vVelocityGradient_x(Mesh,Boundaries)

%% velocity boundary conditions on left and right
[Dvx,bx_bc] = Gradient_Nonuniform_Dirichlet(Mesh.Nx, ...
                                     Mesh.v_dx_centre, ...
                                     Mesh.v_dx,...
                                     Boundaries.Left_v, ...
                                     Boundaries.Right_v);


Iy = speye(Mesh.Ny+1,Mesh.Ny+1);

Dv = kron(Iy,Dvx);

bv_bc = zeros(Mesh.Nx*(Mesh.Ny+1),1);

LeftWall = index(1:Mesh.Nx,1,Mesh.Nx);
RightWall = index(1:Mesh.Nx,Mesh.Ny+1,Mesh.Nx);

for n= 1:length(LeftWall)
   bv_bc(LeftWall(n)) = bx_bc(1);
   bv_bc(RightWall(n)) = bx_bc(Mesh.Nx);
end


end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end
