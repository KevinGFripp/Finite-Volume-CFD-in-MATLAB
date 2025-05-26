function div_v_y = Divergence_y(Mesh)
%% divergence of P-grid face-interpolated v velocity
Dy = Gradient_CellWise_Nonuniform(Mesh.Ny+1,Mesh.dy);

Ix = speye(Mesh.Nx,Mesh.Nx);

div_v_y = kron(Dy,Ix);

end