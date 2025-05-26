function PlotMeshCentres(Mesh)
Nx = Mesh.Nx;
Ny = Mesh.Ny;

ugridx = zeros(Nx+1,Ny);
ugridy = zeros(Nx+1,Ny);
vgridx = zeros(Nx,Ny+1);
vgridy = zeros(Nx,Ny+1);
Pgridx = zeros(Nx,Ny);
Pgridy = zeros(Nx,Ny);

for x = 1 :Nx+1
    for y = 1:Ny
ugridx(x,y) = Mesh.u_centre_x(x);
ugridy(x,y) = Mesh.u_centre_y(y);
    end
end

for x = 1 :Nx
    for y = 1:Ny+1
vgridx(x,y) = Mesh.v_centre_x(x);
vgridy(x,y) = Mesh.v_centre_y(y);
    end
end

for x = 1 :Nx
    for y = 1:Ny
Pgridx(x,y) = Mesh.P_centre_x(x);
Pgridy(x,y) = Mesh.P_centre_y(y);
    end
end

plot(reshape(Pgridx,Nx*Ny,1),reshape(Pgridy,Nx*Ny,1), ...
     'ko','MarkerFaceColor','k',MarkerSize=5./sqrt(Nx));
hold on
plot(reshape(ugridx,(Nx+1)*Ny,1),reshape(ugridy,(Nx+1)*Ny,1), ...
    'b>','MarkerFaceColor','b',Markersize=5/sqrt(Nx));
plot(reshape(vgridx,Nx*(Ny+1),1),reshape(vgridy,Nx*(Ny+1),1), ...
    'r^','MarkerFaceColor','r',MarkerSize=5/sqrt(Nx));
hold off


end