function [x,rnorm,iter] = MGBicgstab_BlockMatrix(Grid,b,x0,maxiter,TOL,BLOCKSIZE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multigrid Preconditioned Bi-Conjugate Gradient Stabilised Method %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% random initial residual improves stability
rhat = rand(length(b),1);

x = x0;
r = b-Grid(1).A*x;
rho = dot(rhat,r);
P = r;
rnorm = zeros(maxiter,1);

y = zeros(length(b),1);
z = zeros(length(b),1);

for I = 1:maxiter
iter = I;

%% Apply left-preconditioning
[y,~,~] = MultiGrid_BlockMatrix_PreAllocatedGrid(Grid, ...
                                                 P, ...
                                                 zeros(length(b),1), ...
                                                 1, ...
                                                 1e-7, ...
                                                 'F', ...
                                                 BLOCKSIZE);
v = Grid(1).A*y;
alpha = rho/dot(rhat,v);
h = x +alpha*y;
s = r -alpha*v;

rnorm(I) = sqrt(sum(s.^2))/length(s);
if (rnorm(I) <= TOL)
    x = h;
    rnorm = squeeze(rnorm);
    break
end

%% Apply left-preconditioning
[z,~,~] = MultiGrid_BlockMatrix_PreAllocatedGrid(Grid, ...
                                                 s, ...
                                                 zeros(length(b),1), ...
                                                 1, ...
                                                 1e-7, ...
                                                 'F', ...
                                                 BLOCKSIZE);

t = Grid(1).A*z;
omega = dot(t,s)/dot(t,t);
x = h+omega*z;
r = s-omega*t;

rnorm(I) = sqrt(sum(r.^2))/length(r);
if (rnorm(I) <= TOL)
    rnorm = squeeze(rnorm);
    break;
end

rho1 = dot(rhat,r);
Beta = (rho1/rho)*(alpha/omega);
P = r +Beta*(P - omega*v);
rho = rho1;

end


end