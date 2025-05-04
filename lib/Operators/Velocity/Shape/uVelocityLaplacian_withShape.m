function [Lu,bu_bc] = uVelocityLaplacian_withShape(Mesh,Boundaries,Shape)

%% velocity defined on left and right boundaries
 Lux = Laplacian_Nonuniform(Mesh.Nx+1,Mesh.u_dx_centre,Mesh.u_dx);


%% velocity boundary conditions top and bottom
[Luy,by_bc] = Laplacian_Nonuniform_Dirichlet(Mesh.Ny, ...
                                     Mesh.u_dy_centre, ...
                                     Mesh.u_dy, ...
                                     Boundaries.Bottom_u, ...
                                     Boundaries.Top_u);

Iy = speye(Mesh.Ny,Mesh.Ny);
Ix = speye(Mesh.Nx+1,Mesh.Nx+1);

% Left and Right
L_left_right =  kron(Iy,Lux);

% Bottom and Top
L_bottom_top = kron(Luy,Ix);

Lu = L_bottom_top + L_left_right;


bu_bc = zeros((Mesh.Nx+1)*Mesh.Ny,1);

BottomWall = index(1:Mesh.Nx+1,1,Mesh.Nx+1);
TopWall = index(1:Mesh.Nx+1,Mesh.Ny,Mesh.Nx+1);

for n= 1:length(BottomWall)
   bu_bc(BottomWall(n)) = by_bc(1);
   bu_bc(TopWall(n)) = by_bc(Mesh.Ny);
end

%% Shape 
igrid = Shape.u_grid_i;
jgrid = Shape.u_grid_j;

%% Shape Interior

for i = 2:(length(igrid)-1)
    for j = 1:length(jgrid)
    ind = index(igrid(i),jgrid(j),Mesh.Nx+1);
    Lu(ind,:) = 0;
    end
end

%% Left wall boundary
% keep y derivative, change x derivative
m = igrid(1);

for j = 1:length(jgrid)
  ind = index(m,jgrid(j),Mesh.Nx+1);

  stencil = [(Mesh.u_dx_centre(m-3)+Mesh.u_dx_centre(m-2) +Mesh.u_dx_centre(m-1)) ...
             (Mesh.u_dx_centre(m-2) +Mesh.u_dx_centre(m-1)) ...
              Mesh.u_dx_centre(m-1)];
  Cm = SolveForCoefficient_RightEdge(stencil,Mesh.u_dx(m));
  coeff = [Cm(1) Cm(2) Cm(3) -sum(Cm)];
  inds = index((m-3):(m),jgrid(j),Mesh.Nx+1);

  Lu(ind,:) = L_bottom_top(ind,:);
  Lu(ind,inds) =Lu(ind,inds) + coeff;
end

%% One-from Left wall boundary
m = igrid(1)-1;
% keep y derivative, change x derivative

for j = 1 :length(jgrid)
  ind = index(m,jgrid(j),Mesh.Nx+1);

  stencil = [(Mesh.u_dx_centre(m-2) +Mesh.u_dx_centre(m-1)) ...
             (Mesh.u_dx_centre(m-1)) ...
              Mesh.u_dx_centre(m)];
  Cm = SolveForCoefficient_OneFromRightEdge(stencil,Mesh.u_dx(m));
  coeff = [Cm(1) Cm(2) -sum(Cm) Cm(3)];
  inds = index((m-2):(m+1),jgrid(j),Mesh.Nx+1);

  Lu(ind,:) = L_bottom_top(ind,:);
  Lu(ind,inds) =Lu(ind,inds) + coeff;
end

%% Right Boundary
% keep y derivative, change x derivative
m = igrid(end);

for j = 1 :length(jgrid)
  ind = index(m,jgrid(j),Mesh.Nx+1);
  
  stencil = [(Mesh.u_dx_centre(m-3)+Mesh.u_dx_centre(m-2) +Mesh.u_dx_centre(m-1)) ...
             (Mesh.u_dx_centre(m-2) +Mesh.u_dx_centre(m-1)) ...
              Mesh.u_dx_centre(m-1)];
  Cm = SolveForCoefficient_RightEdge(stencil,Mesh.u_dx(m));
  coeff = [Cm(1) Cm(2) Cm(3) -sum(Cm)];
  inds = index((m-3):(m),jgrid(j),Mesh.Nx+1);

  Lu(ind,:) = L_bottom_top(ind,:);
  Lu(ind,inds) =Lu(ind,inds) + coeff;
end

%% Bottom Boundary
% keep x derivative, change y derivative

m = jgrid(1)-1;

for i = 1 :length(igrid)
  ind = index(igrid(i),m,Mesh.Nx+1);

  stencil = [(Mesh.u_dy_centre(m-2) +Mesh.u_dy_centre(m-1))...
             Mesh.u_dy_centre(m-1) Mesh.u_dy(m)];
  Cm = SolveForCoefficient_OneFromRightEdge(stencil,Mesh.u_dy(m));
  coeff = [Cm(1) Cm(2) -(sum(Cm)+Cm(3))];
  inds = index(igrid(i),(m-2):(m),Mesh.Nx+1);

  Lu(ind,:) = L_left_right(ind,:);
  Lu(ind,inds) =Lu(ind,inds) + coeff;
end

%% One from bottom boundary
% keep x derivative, change y derivative

m = jgrid(1)-2;

for i = 1 :length(igrid)
  ind = index(igrid(i),m,Mesh.Nx+1);

  stencil = [(Mesh.u_dy_centre(m-2) +Mesh.u_dy_centre(m-1)) ...
            (Mesh.u_dy_centre(m-1)) ...
            Mesh.u_dy_centre(m)];
  Cm = SolveForCoefficient_OneFromRightEdge(stencil,Mesh.u_dy(m));
  coeff = [Cm(1) Cm(2) -sum(Cm) Cm(3)];
  inds = index(igrid(i),(m-2):(m+1),Mesh.Nx+1);
 
  Lu(ind,:) = L_left_right(ind,:);
  Lu(ind,inds) =Lu(ind,inds) + coeff;
end

%% Top boundary
% keep x derivative, change y derivative
m = jgrid(end)+1;

for i = 1 :length(igrid)
  ind = index(igrid(i),m,Mesh.Nx+1);

  stencil = [Mesh.u_dy(m) Mesh.u_dy_centre(m) ...
            (Mesh.u_dy_centre(m) +Mesh.u_dy_centre(m+1))];
  Cm = SolveForCoefficient(stencil,Mesh.u_dy(m)); 
  coeff = [-(sum(Cm)+Cm(1)) Cm(2) Cm(3)];
  inds = index(igrid(i),(m):(m+2),Mesh.Nx+1);

  Lu(ind,:) = L_left_right(ind,:);
  Lu(ind,inds) =Lu(ind,inds) + coeff;
end


end

function C = SolveForCoefficient(stencil,dx)
% %% 4 point stencil, 2nd order accurate
% % -ah0 +bh1 +ch2 = 0
% % -ah0^3 +bh1^3 +ch2^3 = 0
% % +ah0^2 +bh1^2 +ch2^2 = 2h0^2
CoefficientMatrix = zeros(3,3);
 b = zeros(3,1);
 b(3,1) = 2 * dx.^2;

 CoefficientMatrix(1,:) =[-stencil(1) stencil(2) stencil(3)];
 CoefficientMatrix(2,:) =[-stencil(1) stencil(2) stencil(3)].^3;
 CoefficientMatrix(3,:) =[stencil(1) stencil(2) stencil(3)].^2;
 
  C = CoefficientMatrix\b;
  C = C./dx./dx;
end

function C = SolveForCoefficient_OneFromRightEdge(stencil,dx)
% %% 4 point stencil, 2rd order accurate
% % -ah0 -bh1 +ch2 = 0
% % -ah0^3 -bh1^3 +ch2^3 = 0
% % +ah0^2 +bh1^2 +ch2^2 = 2h0^2
CoefficientMatrix = zeros(3,3);
 b = zeros(3,1);
 b(3,1) = 2 * dx.^2;

 CoefficientMatrix(1,:) =[-stencil(1) -stencil(2) stencil(3)];
 CoefficientMatrix(2,:) =[-stencil(1) -stencil(2) stencil(3)].^3;
 CoefficientMatrix(3,:) =[stencil(1) stencil(2) stencil(3)].^2;

  C = CoefficientMatrix\b;
  C = C./dx./dx;
end

function C = SolveForCoefficient_RightEdge(stencil,dx)
% %% 4 point stencil, 2rd order accurate
% % -ah0 -bh1 +ch2 = 0
% % -ah0^3 -bh1^3 +ch2^3 = 0
% % +ah0^2 +bh1^2 +ch2^2 = 2h0^2
CoefficientMatrix = zeros(3,3);
 b = zeros(3,1);
 b(3,1) = 2 * dx.^2;

 CoefficientMatrix(1,:) =[-stencil(1) -stencil(2) -stencil(3)];
 CoefficientMatrix(2,:) =[-stencil(1) -stencil(2) -stencil(3)].^3;
 CoefficientMatrix(3,:) =[stencil(1) stencil(2) stencil(3)].^2;

  C = CoefficientMatrix\b;
  C = C./dx./dx;
end

function k = index(i,j,Nx)

k = i + (j-1)*Nx;

end
