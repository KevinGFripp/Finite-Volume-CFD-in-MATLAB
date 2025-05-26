function Iv = Interp_v_velocity_to_CellCentre(Mesh)

%% forward difference, no boundary conditions
Nx = Mesh.Nx;
Ny = Mesh.Ny;
kernel_y = spdiags([ones(Ny+1,1) ones(Ny+1,1)],[0 1],Ny,Ny+1);
Ix = speye(Nx,Nx);

Iv = kron(kernel_y,Ix)./2;

end