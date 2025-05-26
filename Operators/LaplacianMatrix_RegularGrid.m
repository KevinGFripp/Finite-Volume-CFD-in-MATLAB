function [L,b] = LaplacianMatrix_RegularGrid(Nx,Ny,dx,dy,top)

% Laplacian solved on interior with Boundary Conditions applied to matrix

% Dirichlet boundary conditions on all sides
L = MakeLaplacianMatrix(Nx,Ny,dx,dy);

% Apply Neumann Boundary Conditions on 3 sides
RightWall = index(Nx,1:Ny,Nx);
LeftWall =  index(1,1:Ny,Nx);
BottomWall = index(1:Nx,1,Nx);
TopWall = index(1:Nx,Ny,Nx);

b = zeros(Nx*Ny,1);
b(TopWall) = 2*top/dy/dy;

%% walls not the same length top/bottom versus left/right
for n=1:length(LeftWall)
L(RightWall(n),RightWall(n)) = L(RightWall(n),RightWall(n)) - 1./dx./dx;
L(LeftWall(n),LeftWall(n)) = L(LeftWall(n),LeftWall(n)) - 1./dx./dx;
end

for n=1:length(TopWall)
L(BottomWall(n),BottomWall(n)) = L(BottomWall(n),BottomWall(n)) - 1./dy./dy;
L(TopWall(n),TopWall(n)) = L(TopWall(n),TopWall(n)) - 1./dy./dy;
end


%% ----------------------------


end

function k = index(i,j,Nx)

k = i +(j-1)*Nx;

end


function L = MakeLaplacianMatrix(Nx,Ny,dx,dy)

Ln = spdiags([ones(Nx,1) -2.*ones(Nx,1) ones(Nx,1)],-1:1,Nx,Nx)./dx./dx;
Lm = spdiags([ones(Ny,1) -2.*ones(Ny,1) ones(Ny,1)],-1:1,Ny,Ny)./dy./dy;
Identitym =speye(Ny);
Identityn =speye(Nx);
L = kron(Lm,Identityn)+kron(Identitym,Ln);

end