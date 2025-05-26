function Iu = Interp_u_velocity_to_CellCentre(Mesh)

%% forward interpolation, no boundary conditions
Nx = Mesh.Nx;
Ny = Mesh.Ny;

kernel_x = spdiags([ones(Nx+1,1) ones(Nx+1,1)],[0 1],Nx,Nx+1);
Iy = speye(Ny,Ny);

Iu = kron(Iy,kernel_x)./2;

end