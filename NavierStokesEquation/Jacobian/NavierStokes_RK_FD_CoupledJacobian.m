function [J_FD] = NavierStokes_RK_FD_CoupledJacobian(uvector,vvector, ...
                                            uOperator,vOperator, ...
                                            Mesh,Boundaries,Parameters)


%% Finite difference Jacobian approximation
 Nx = Mesh.Nx;
 Ny = Mesh.Ny;
 JCOLS = (Nx+1)*Ny + Nx*(Ny+1);
 eps = 1e-7;
 J_FD = zeros(JCOLS,JCOLS);
 eps_u_array = zeros((Nx+1)*Ny,1);
 eps_v_array = zeros(Nx*(Ny+1),1);

 [du_dt,dv_dt] = NavierStokes_RK_NoP(uvector,vvector, ...
                                 uOperator,vOperator, ...
                                 Mesh,Boundaries,Parameters);

%% uvvector
for n = 1:(Nx+1)*Ny
    temp = eps_u_array;
    temp(n) =eps;
    [du_dt_eps,dv_dt_eps] = NavierStokes_RK_NoP(uvector+temp,vvector, ...
                                            uOperator,vOperator, ...
                                            Mesh,Boundaries,Parameters);

    J_FD(:,n) = ([du_dt_eps;dv_dt_eps] - [du_dt;dv_dt])./eps;
end

%% vvvector
for n = 1:Nx*(Ny+1)
    temp = eps_v_array;
    temp(n) =eps;
    [du_dt_eps,dv_dt_eps] = NavierStokes_RK_NoP(uvector,vvector+temp, ...
                                                uOperator,vOperator, ...
                                                Mesh,Boundaries,Parameters);

    J_FD(:,n+(Nx+1)*Ny) = ([du_dt_eps;dv_dt_eps] - [du_dt;dv_dt])./eps;
end

J_FD = sparse(J_FD);
   
end
