function div_v_y = Divergence_y_withShape(Mesh,Shape)
%% divergence of P-grid face-interpolated v velocity
% Equal to 0 on the interior of the shape
Dy = Gradient_CellWise_Nonuniform(Mesh.Ny+1,Mesh.dy);

Ix = speye(Mesh.Nx,Mesh.Nx);

div_v_y = kron(Dy,Ix);

%% Shape interior
igrid = Shape.p_grid_i;
jgrid = Shape.p_grid_j;

for i=1:length(igrid)
    for j=1:length(jgrid)
    ind = index(igrid(i),jgrid(j),Mesh.Nx);
    div_v_y(ind,:) = 0;
    end
end

end