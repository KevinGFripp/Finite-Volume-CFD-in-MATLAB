function Grid = PreallocateHierachyOfGrids(A, ...
                                           b, ...
                                           x0, ...
                                           GRIDS, ...
                                           Nx, ...
                                           Ny, ...
                                           BLOCKSIZE)

%%A matrix on grid n
%%x solution on grid n
%%b vector on grid 
%%R restriction operator from grid n to grid n+1
%%I interpolation operator from grid n+1 to n
%%Nx,Ny size of grid n

%%Initialise structs
for n=1:GRIDS
Grid(n) =struct('A',1,'x',1,'b',1,'R',1,'I',1,'Nx',1,'Ny',1);
end

%% Fine mesh
Grid(1) =struct('A',A,'x',x0,'b',b,'R',1,'I',1,'Nx',Nx,'Ny',Ny);

[Grid(1).R,...
 Grid(1).I,...
 Grid(2).Nx,...
 Grid(2).Ny] = RestrictionandProlongationOperator_BM(Nx,Ny,BLOCKSIZE);

%% Coarser Grids
for n=2:GRIDS  

 if (n < GRIDS)
[Grid(n).R,...
 Grid(n).I,...
 Grid(n+1).Nx,...
 Grid(n+1).Ny] = RestrictionandProlongationOperator_BM(Grid(n).Nx,Grid(n).Ny,BLOCKSIZE);
 end

Grid(n).A = Grid(n-1).R*(Grid(n-1).A)*Grid(n-1).I;
Grid(n).x = zeros(BLOCKSIZE*(Grid(n).Nx)*(Grid(n).Ny),1);

end

end

function [R_2D,I_2D,Nx_2h,Ny_2h] = RestrictionandProlongationOperator_BM(Nx,Ny,BLOCKSIZE)

%%Odd or even kernels
K_N_odd = [1,2,1];
K_N_even = [1,3,3,1];

Nx_2h = (ceil(Nx/2)-1);
Ny_2h = (ceil(Ny/2)-1);

Rx_2h =zeros(Nx_2h,Nx);
Ry_2h =zeros(Ny_2h,Ny);


%%Block kernel
kernel_block =diag(ones(BLOCKSIZE,1));

%% Choose stencil based upon whether odd or even N
 if(mod(Nx,2) == 0)
  stencilx = K_N_even;
  R_2h_factorx = 1/8;
 else
  stencilx = K_N_odd;
  R_2h_factorx = 1/4;
 end

 if(mod(Ny,2) == 0)
  stencily = K_N_even;
  R_2h_factory = 1/8;
 else
  stencily = K_N_odd;
  R_2h_factory = 1/4;
end

%%Make restriction matrix
for n=1:Nx_2h
Rx_2h(n,(1+2*(n-1)):(2*(n-1)+length(stencilx)))=stencilx;
end
for n=1:Ny_2h
Ry_2h(n,(1+2*(n-1)):(2*(n-1)+length(stencily)))=stencily;
end

Rx_2h = sparse(Rx_2h);
Ry_2h = sparse(Ry_2h);

%%2D restriction is Kronecker product kron(R_2h,R_2h)
  R_2D = kron(R_2h_factorx.*Rx_2h,R_2h_factory.*Ry_2h);
  R_2D = kron(R_2D,kernel_block);
  I_2D = transpose(R_2D);


end