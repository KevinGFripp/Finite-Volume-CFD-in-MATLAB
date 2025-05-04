function div_u_x = Divergence_x(Mesh)
%% divergence of P-grid face-interpolated u velocity
Dx = Gradient_CellWise_Nonuniform(Mesh.Nx+1,Mesh.dx);

Iy = speye(Mesh.Ny,Mesh.Ny);

div_u_x = kron(Iy,Dx);

end