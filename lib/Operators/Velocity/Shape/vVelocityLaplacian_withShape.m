function [Lv,bv_bc] = vVelocityLaplacian_withShape(Mesh,Boundaries,Shape)


%% velocity boundary conditions on left and right
[Lvx,bx_bc] = Laplacian_Nonuniform_Dirichlet(Mesh.Nx, ...
                                     Mesh.v_dx_centre, ...
                                     Mesh.v_dx,...
                                     Boundaries.Left_v, ...
                                     Boundaries.Right_v);

[Lvx_Nm] = Laplacian_Nonuniform_Neumann(Mesh.Nx, ...
                                     Mesh.v_dx_centre, ...
                                     Mesh.v_dx);

if(strcmp(Boundaries.RightWallType,'Free'))
  Lvx(Mesh.Nx,:) = Lvx_Nm(Mesh.Nx,:);
end

%% velocity defined on top and bottom
Lvy = Laplacian_Nonuniform(Mesh.Ny+1, ...
                           Mesh.v_dy_centre, ...
                           Mesh.v_dy);

Iy = speye(Mesh.Ny+1,Mesh.Ny+1);
Ix = speye(Mesh.Nx,Mesh.Nx);

% Neumann left-right
L_left_right_Nm = kron(Iy,Lvx_Nm);

% Left and Right
L_left_right =  kron(Iy,Lvx);

% Bottom and Top
L_bottom_top = kron(Lvy,Ix);

Lv = L_left_right + L_bottom_top;

bv_bc = zeros(Mesh.Nx*(Mesh.Ny+1),1);

LeftWall = index(1,1:Mesh.Ny+1,Mesh.Nx);
RightWall = index(Mesh.Nx,1:Mesh.Ny+1,Mesh.Nx);

for n= 1:length(LeftWall)
   bv_bc(LeftWall(n)) = bx_bc(1);
   bv_bc(RightWall(n)) = bx_bc(Mesh.Nx);
end



%% Shape 
igrid = Shape.v_grid_i;
jgrid = Shape.v_grid_j;

%% Shape Interior

for i = 1 : length(igrid)
    for j = 2 : (length(jgrid) - 1)
    ind = index(igrid(i),jgrid(j),Mesh.Nx);
    Lv(ind,:) = 0;
    end
end

%% Left wall boundary
% keep y derivative, change x derivative
m = igrid(1)-1;

for j = 1:length(jgrid)
  ind = index(m,jgrid(j),Mesh.Nx);

  stencil = [(Mesh.v_dx_centre(m-2) +Mesh.v_dx_centre(m-1))...
             Mesh.v_dx_centre(m-1) Mesh.v_dx(m)];
  Cm = SolveForCoefficient_OneFromRightEdge(stencil,Mesh.v_dx(m));
  coeff = [Cm(1) Cm(2) -(sum(Cm)+Cm(3))];
  inds = index((m-2):(m),jgrid(j),Mesh.Nx);
  
  Lv(ind,:) = L_bottom_top(ind,:);
  Lv(ind,inds) =Lv(ind,inds) + coeff;
end

%% One-from Left wall boundary
m = igrid(1)-2;
% keep y derivative, change x derivative

for j = 1:length(jgrid)
  ind = index(m,jgrid(j),Mesh.Nx);

  stencil = [(Mesh.v_dx_centre(m-2) +Mesh.v_dx_centre(m-1)) ...
             (Mesh.v_dx_centre(m-1)) ...
              Mesh.v_dx_centre(m)];
  Cm = SolveForCoefficient_OneFromRightEdge(stencil,Mesh.v_dx(m));
  coeff = [Cm(1) Cm(2) -sum(Cm) Cm(3)];
  inds = index((m-2):(m+1),jgrid(j),Mesh.Nx);

  Lv(ind,:) = L_bottom_top(ind,:);
  Lv(ind,inds) =Lv(ind,inds) + coeff;
end

%% Right Boundary
% keep y derivative, change x derivative
m = igrid(end)+1;

for j = 1:length(jgrid)
  ind = index(m,jgrid(j),Mesh.Nx);

  stencil = [Mesh.v_dx(m) Mesh.v_dx_centre(m) ...
            (Mesh.v_dx_centre(m) +Mesh.v_dx_centre(m+1))];
  Cm = SolveForCoefficient(stencil,Mesh.v_dx(m)); 
  coeff = [-(sum(Cm)+Cm(1)) Cm(2) Cm(3)];
  inds = index((m):(m+2),jgrid(j),Mesh.Nx);

  Lv(ind,:) = L_bottom_top(ind,:);
  Lv(ind,inds) =Lv(ind,inds) + coeff;
end

%% Bottom Boundary
% keep x derivative, change y derivative
m = jgrid(1);

for i = 1 :length(igrid)
  ind = index(igrid(i),m,Mesh.Nx);

  stencil = [(Mesh.v_dy_centre(m-3)+Mesh.v_dy_centre(m-2) +Mesh.v_dy_centre(m-1)) ...
             (Mesh.v_dy_centre(m-2) +Mesh.v_dy_centre(m-1)) ...
              Mesh.v_dy_centre(m-1)];
  Cm = SolveForCoefficient_RightEdge(stencil,Mesh.v_dy(m));
  coeff = [Cm(1) Cm(2) Cm(3) -sum(Cm)];
  inds = index(igrid(i),(m-3):(m),Mesh.Nx);
  
  Lv(ind,:) = L_left_right(ind,:);
  Lv(ind,inds) =Lv(ind,inds) + coeff;
end

%% One from bottom boundary
% keep x derivative, change y derivative
m = jgrid(1)-1;

for i = 1 :length(igrid)
  ind = index(igrid(i),m,Mesh.Nx);

  stencil = [(Mesh.v_dy_centre(m-2) +Mesh.v_dy_centre(m-1)) ...
             (Mesh.v_dy_centre(m-1)) ...
              Mesh.v_dy_centre(m)];
  Cm = SolveForCoefficient_OneFromRightEdge(stencil,Mesh.v_dy(m));
  coeff = [Cm(1) Cm(2) -sum(Cm) Cm(3)];
  inds = index(igrid(i),(m-2):(m+1),Mesh.Nx);

  Lv(ind,:) = L_left_right(ind,:);
  Lv(ind,inds) =Lv(ind,inds) + coeff;
end

%% Top boundary
% keep x derivative, change y derivative
m = jgrid(end);

for i = 1 :length(igrid)
  ind = index(igrid(i),m,Mesh.Nx);

  stencil = [(Mesh.v_dy_centre(m-3)+Mesh.v_dy_centre(m-2) +Mesh.v_dy_centre(m-1)) ...
             (Mesh.v_dy_centre(m-2) +Mesh.v_dy_centre(m-1)) ...
              Mesh.v_dy_centre(m-1)];
  Cm = SolveForCoefficient_RightEdge(stencil,Mesh.v_dy(m));
  coeff = [Cm(1) Cm(2) Cm(3) -sum(Cm)];
  inds = index(igrid(i),(m-3):(m),Mesh.Nx);

  Lv(ind,:) = L_left_right(ind,:);
  Lv(ind,inds) =Lv(ind,inds) + coeff;
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
