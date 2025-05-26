function div_u_x = Divergence_x_withShape(Mesh,Shape)
%% divergence of P-grid face-interpolated u velocity
% Equal to 0 on the interior of the shape

Dx = Gradient_CellWise_Nonuniform(Mesh.Nx+1,Mesh.dx);

Iy = speye(Mesh.Ny,Mesh.Ny);

div_u_x = kron(Iy,Dx);

%% Shape interior
igrid = Shape.p_grid_i;
jgrid = Shape.p_grid_j;

for i=1:length(igrid)
    for j=1:length(jgrid)
    ind = index(igrid(i),jgrid(j),Mesh.Nx);
    div_u_x(ind,:) = 0;
    end
end


end

function k = index(i,j,Nx)

k = i +(j-1)*Nx;

end
