function Mesh = MakeRectilinearMesh(Nx,Ny,Width,Height,growthx,growthy,Shapes)

Mesh.Nx = Nx;
Mesh.Ny = Ny;
Mesh.W = Width;
Mesh.H = Height;

dx = Width/Nx;
dy = Height/Ny;

% Pressure Grid and faces
Mesh.P_centre_x = zeros(Nx,1);
Mesh.P_centre_y = zeros(Ny,1);
Mesh.P_dx_centre = zeros(Nx-1,1);
Mesh.P_dy_centre = zeros(Ny-1,1);
Mesh.dx = zeros(Nx,1);
Mesh.dy = zeros(Ny,1);
Mesh.xface = zeros(Nx+1,1);
Mesh.yface = zeros(Ny+1,1);


% velocity cell sizes -> Distance between P centroids
Mesh.u_centre_x = zeros(Nx+1,1);
Mesh.u_centre_y = zeros(Ny,1);
Mesh.v_centre_x = zeros(Nx,1);
Mesh.v_centre_y = zeros(Ny+1,1);
Mesh.u_dx = zeros(Nx+1,1);
Mesh.u_dy = zeros(Ny,1);
Mesh.v_dx = zeros(Nx,1);
Mesh.v_dy = zeros(Ny+1,1);
Mesh.u_dx_centre = zeros(Nx+1-1,1);
Mesh.u_dy_centre = zeros(Ny-1,1);
Mesh.v_dx_centre = zeros(Nx-1,1);
Mesh.v_dy_centre = zeros(Ny+1-1,1);




%% Define series of boundary intervals
if(isempty(Shapes))
    NUMSHAPES = 0;
else
NUMSHAPES = length(Shapes);
end
Intervalsx = zeros(2*NUMSHAPES+2,1);
Intervalsy = zeros(2*NUMSHAPES+2,1);

 for S = 1 : NUMSHAPES
 
   Intervalsx(S+1,1) = Shapes(S).p_grid_i(1)+1;
   Intervalsx(S+2,1) = Shapes(S).p_grid_i(end)+1;
 
   Intervalsy(S+1,1) = Shapes(S).p_grid_j(1)+1;
   Intervalsy(S+2,1) = Shapes(S).p_grid_j(end)+1;
 
 end
 
 Intervalsx(1) = 1;
 Intervalsx(2*NUMSHAPES+2) = Nx+1;
 Intervalsy(1) = 1;
 Intervalsy(2*NUMSHAPES+2) = Ny+1;


%% Construct mesh regions
for m = 1:(2*NUMSHAPES+1)

 x0 = dx*(Intervalsx(m)-1);
 y0 = dy*(Intervalsy(m)-1);
 indm_x = Intervalsx(m):Intervalsx(m+1);
 indm_y = Intervalsy(m):Intervalsy(m+1);

  Nm_x = length(indm_x)-1;
  Nm_y = length(indm_y)-1;

  Width_n = Nm_x * dx;
  Height_n = Nm_y * dy;

  Mesh.xface(indm_x) = x0 + CreateNonUniformMesh_x(Width_n,Nm_x,growthx(m));
  Mesh.yface(indm_y) = y0 + CreateNonUniformMesh_y(Height_n,Nm_y,growthy(m));

end


Mesh.dx = Mesh.xface(2:Nx+1)-Mesh.xface(1:Nx);
Mesh.dy = Mesh.yface(2:Ny+1)-Mesh.yface(1:Ny);

% Pressure cell centroids
Mesh.P_centre_x = 0.5*(Mesh.xface(2:Nx+1) + Mesh.xface(1:Nx));
Mesh.P_centre_y = 0.5*(Mesh.yface(2:Ny+1) + Mesh.yface(1:Ny));

% cell distance between cell centres
Mesh.P_dx_centre = Mesh.P_centre_x(2:Nx) - Mesh.P_centre_x(1:Nx-1);
Mesh.P_dy_centre = Mesh.P_centre_y(2:Ny) - Mesh.P_centre_y(1:Ny-1);

% velocity non-staggered cell sizes
Mesh.v_dx = Mesh.dx;
Mesh.u_dy = Mesh.dy;

%  u velocity staggered direction cell size
Mesh.u_dx(2:Nx) = Mesh.P_centre_x(2:Nx) - Mesh.P_centre_x(1:Nx-1);
% boundaries adopt Pressure grid cell size
Mesh.u_dx(1) = Mesh.dx(1);
Mesh.u_dx(Nx+1) = Mesh.dx(Nx);

% velocity cell centre
Mesh.u_centre_y = Mesh.P_centre_y;
Mesh.u_centre_x(1) = 0;
Mesh.u_centre_x(Nx+1) = Width;
Mesh.u_centre_x(2:Nx) = 0.5*(Mesh.P_centre_x(2:Nx) + Mesh.P_centre_x(1:Nx-1));

% cell distance between cell centres
Mesh.u_dx_centre = Mesh.u_centre_x(2:Nx+1) - Mesh.u_centre_x(1:Nx);
Mesh.u_dy_centre = Mesh.u_centre_y(2:Ny) - Mesh.u_centre_y(1:Ny-1);

%  v velocity staggered direction cell size
Mesh.v_dy(2:Ny) = Mesh.P_centre_y(2:Ny) - Mesh.P_centre_y(1:Ny-1);
% boundaries adopt Pressure grid cell size
Mesh.v_dy(1) = Mesh.dy(1);
Mesh.v_dy(Ny+1) = Mesh.dy(Ny);

% velocity cell centre
Mesh.v_centre_x = Mesh.P_centre_x;
Mesh.v_centre_y(1) = 0;
Mesh.v_centre_y(Ny+1) = Height;
Mesh.v_centre_y(2:Ny) = 0.5*(Mesh.P_centre_y(2:Ny) + Mesh.P_centre_y(1:Ny-1));

% cell distance between cell centres
Mesh.v_dx_centre = Mesh.v_centre_x(2:Nx) - Mesh.v_centre_x(1:Nx-1);
Mesh.v_dy_centre = Mesh.v_centre_y(2:Ny+1) - Mesh.v_centre_y(1:Ny);


end

function xn = CreateNonUniformMesh_x(Width_x,Nn,growthx)

indx = 0:Nn;
dxi = 2/Nn;
xi = (dxi * indx) - 1;
xn = GrowthRate(growthx,xi,Width_x);

end

function yn = CreateNonUniformMesh_y(Height_n,Nn,growthy)

indy = 0:Nn;
dyi = 2/Nn;
yi = (dyi * indy) - 1;
yn = GrowthRate(growthy,yi,Height_n);

end

function xnew = GrowthRate(wx,xi,W)

xnew = (W/2) * (1 + tanh(wx*xi)./tanh(wx));

end