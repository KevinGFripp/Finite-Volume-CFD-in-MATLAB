function [x,rnorm,iter] = MultiGrid_BlockMatrix_PreAllocatedGrid(Grid, ...
                                                                 b, ...
                                                                 x0, ...
                                                                 maxiterations, ...
                                                                 TOL, ...
                                                                 type, ...
                                                                 BLOCKSIZE)
%%%%%%%%%%%%%%%%%%%%
%%%% Multigrid %%%%%
%%%% F | V cycle %%%
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
%%For block matrices%%
%%%%%%%%%%%%%%%%%%%%%%

%% A = sparse square block matrix from Nx x Ny grid of Laplacian type 
%% Blocks of size BLOCKSIZE x BLOCKSIZE
%% b = Right-Hand side solution vector
%% x0 = initial solution
%% TOL = error tolerance
%% type = 'V','F' cycle

%%%%%%Setup%%%%%%%%
%%%%%%%%%%%%%%%%%%%

if(isempty(TOL))
  TOL = 1e-6;
end

if(isempty(maxiterations))
    maxiterations = 10;
end

GRIDS = length(Grid);

%% Choose Cycle
if(strcmp(type,'V'))
  F_Cycles = 0;
else
  F_Cycles = GRIDS-2;
end

%% Update RHS and initial solution
Grid(1).x = x0;
Grid(1).b = b;

%% Error
rnorm=zeros(maxiterations,1);

%% Iterations required for TOL
iter = 0;

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%



for I = 1 : maxiterations
    
   %% fine grid
   %% pre-smooth
   [Grid(1).x,~] = GaussSeidel(Grid(1).A, ...
                               Grid(1).b, ...
                               Grid(1).x, ...
                               12);

   %% restrict to coarser mesh
   Grid(2).b = Restrict(Grid(1).R,Residual(Grid(1).A, ...
                                           Grid(1).x, ...
                                           Grid(1).b));
    %% First V Cycle
     for n = 2:(GRIDS-1)
       %% smooth error equation A * xn = r
       [Grid(n).x,~] = GaussSeidel(Grid(n).A, ...
                                   Grid(n).b, ...                               
                                   zeros(BLOCKSIZE*(Grid(n).Nx)*(Grid(n).Ny),1), ...
                                   10);
       %% restrict to coarser mesh
        Grid(n+1).b = Restrict(Grid(n).R,Residual(Grid(n).A, ...
                                                  Grid(n).x, ...
                                                  Grid(n).b));
     end

    %% coarsest grid
    %% solve on coarsest grid
    Grid(GRIDS).x = (Grid(GRIDS).A)\Grid(GRIDS).b;

   
    %% F Cycle loop
    for F=1:F_Cycles
    %% start on lowest grid, interpolate from lowest grid GRIDS to grid GRIDS-F
       for n = (GRIDS-1):-1:(GRIDS-F)

          %% interpolate solution
          Grid(n).x =  Grid(n).x + Prolongate(Grid(n).I,Grid(n+1).x);
          %% post-smooth
          [Grid(n).x,~] = GaussSeidel(Grid(n).A, ...
                                      Grid(n).b, ...
                                      Grid(n).x, ...
                                      5);
       end

       %% Restrict back from grid GRIDS-F to GRIDS
       for n = (GRIDS-F):(GRIDS-1)
       %% smooth error equation A * xn = r
       [Grid(n).x,~] = GaussSeidel(Grid(n).A, ...
                                   Grid(n).b, ...                               
                                   zeros(BLOCKSIZE*(Grid(n).Nx)*(Grid(n).Ny),1), ...
                                   8);
       %% restrict to coarser mesh
        Grid(n+1).b = Restrict(Grid(n).R,Residual(Grid(n).A, ...
                                                  Grid(n).x, ...
                                                  Grid(n).b));
       end

       %% solve on coarsest grid
     Grid(GRIDS).x = (Grid(GRIDS).A)\Grid(GRIDS).b;
 
    end
    %% end of interior F-Cycle
   

     for n = (GRIDS-1):-1:1
       %% interpolate solution
       Grid(n).x =  Grid(n).x + Prolongate(Grid(n).I,Grid(n+1).x);
       %% post-smooth
       [Grid(n).x,~] = GaussSeidel(Grid(n).A, ...
                                   Grid(n).b, ...
                                   Grid(n).x, ...
                                   3);
     end
    
    iter = I;

    %% Error
     r=Residual(Grid(1).A,Grid(1).x,Grid(1).b);
%     rnorm(I) = sqrt(sum(r.^2))./(Grid(1).Nx*Grid(1).Ny);
      rnorm(I) = norm(r)/norm(b);
    if(rnorm(I) <= TOL)        
        break;
    end

end

 x = Grid(1).x;

end

function r = Restrict(R,x)
r = R*x;
end

function p = Prolongate(I,x)
p = I*x;
end

function [x,r] = SOR(A,b,x0,ITER)
%% SOR implementation
w = 0.4;
L = tril(A,-1);
U = triu(A,1);
D = A - L - U;

x = x0;

for n=1:ITER
x = (D + w*L)\(w*b - (w*U +(w-1)*D)*x);
end

r = sqrt(sum((A*x -b).^2));
end

function [x,r] = GaussSeidel(A,b,x0,ITER)
%% forward Gauss-Seidal for pre-smoothing and restriction
L = tril(A,0);
U = triu(A,1);

x = x0;

for n=1:ITER
x = L\(b - U*x);
end

r = sqrt(sum((A*x -b).^2));
end

function [x,r] = BackwardGaussSeidel(A,b,x0,ITER)
%% backward Gauss-Seidel for post-smoothing and intepolation
L = tril(A,-1);
U = triu(A,0);

x = x0;

for n=1:ITER
x = U\(b - L*x);
end

r = sqrt(sum((A*x -b).^2));
end

function r = Residual(A,x,b)
r = b - A*x;
end

