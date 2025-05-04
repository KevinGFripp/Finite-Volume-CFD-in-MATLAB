function [Du,bu_bc] = uVelocityGradient_y(Mesh,Boundaries)

%% velocity boundary conditions on bottom and top
[Duy,by_bc] = Gradient_Nonuniform_Dirichlet(Mesh.Ny, ...
                                     Mesh.u_dy_centre, ...
                                     Mesh.u_dy,...
                                     Boundaries.Bottom_u, ...
                                     Boundaries.Top_u);


Ix = speye(Mesh.Nx+1,Mesh.Nx+1);

Du = kron(Duy,Ix);

bu_bc = zeros((Mesh.Nx+1)*Mesh.Ny,1);

BottomWall = index(1,1:Mesh.Ny,Mesh.Nx+1);
TopWall = index(Mesh.Nx+1,1:Mesh.Ny,Mesh.Nx+1);

for n= 1:length(BottomWall)
   bu_bc(BottomWall(n)) = by_bc(1);
   bu_bc(TopWall(n)) = by_bc(Mesh.Ny);
end


end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end
