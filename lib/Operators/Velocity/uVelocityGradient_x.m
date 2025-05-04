function Du = uVelocityGradient_x(Mesh,Boundaries)



%% velocity doesn't require boundaries
Dux = Gradient_Nonuniform(Mesh.Nx+1, ...
                                  Mesh.u_dx_centre, ...
                                  Mesh.u_dx);


Iy = speye(Mesh.Ny,Mesh.Ny);

Du = kron(Iy,Dux);


end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end
