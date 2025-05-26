function [Solver,uvector_n,vvector_n] = ...
                       NewtonStageSolve(stage,uvector,vvector,...
                                        L,U,Juv,Solver, ...
                                        POperator, ...
                                        uOperator, ...
                                        vOperator, ...
                                        Mesh,Boundaries, ...
                                        Parameters,dt)
%% coefficients
a_s = Solver.a_s;
ci = Solver.ci;
alpha = Solver.a_s(2,1);


%% Stage s
h = ci(stage)*dt;

Un = 0;
Vn = 0;

for k = 1:(stage-1) 
    a_s(stage,k);
    Un = Un + a_s(stage,k)*Solver.du_n_dt(:,k); 
    Vn = Vn + a_s(stage,k)*Solver.dv_n_dt(:,k); 
end

% initial guess
uvector_n = uvector + dt*Un;
vvector_n = vvector + dt*Vn;

[uvector_n,vvector_n] = ApplyVelocityWallBoundaries(uvector_n,vvector_n, ...
                                                    Boundaries,Mesh.Nx,Mesh.Ny);

for n = 1:35
%% Evaluate du/dt
[Solver.du_n_dt(:,stage),...
 Solver.dv_n_dt(:,stage)] = NavierStokes_RK_NoP(uvector_n,vvector_n, ...
                                            uOperator,vOperator, ...
                                            Mesh,Boundaries,Parameters);

%% make F0
F0_u = uvector_n - uvector - dt*Un - dt*alpha*Solver.du_n_dt(:,stage);
F0_v = vvector_n - vvector - dt*Vn - dt*alpha*Solver.dv_n_dt(:,stage);

[F0_u,F0_v] = ApplyVelocityDerivativeWallBoundaries(F0_u,F0_v,...
                                                    Boundaries, ...
                                                    Mesh.Nx,Mesh.Ny);
F0_uv = [F0_u;F0_v];

[deltauv,~,~,~,~] = bicgstabl(Juv,-F0_uv,1e-10,20,L,U);

uvector_n = uvector_n + deltauv(1:(Mesh.Nx+1)*Mesh.Ny);
vvector_n = vvector_n + deltauv(((Mesh.Nx+1)*Mesh.Ny+1):end);

[uvector_n,vvector_n] = ApplyVelocityWallBoundaries(uvector_n,vvector_n, ...
                                                    Boundaries,Mesh.Nx,Mesh.Ny);

res = max(abs(F0_uv));

if(res <= 1e-6) 
    break;
end

end

%% Corrector
[uvector_n,vvector_n,~] = CorrectorStep_RK(uvector_n,vvector_n, ...
                                            POperator,uOperator,vOperator, ...
                                            Mesh,Boundaries,Parameters,h);

[Solver.du_n_dt(:,stage),...
 Solver.dv_n_dt(:,stage)] = NavierStokes_RK_NoP(uvector_n, ...
                                            vvector_n, ...
                                            uOperator,vOperator, ...
                                            Mesh,Boundaries,Parameters);

[Solver.du_n_dt(:,stage),...
 Solver.dv_n_dt(:,stage)] = ApplyVelocityDerivativeWallBoundaries( ...
                              Solver.du_n_dt(:,stage), ...
                              Solver.dv_n_dt(:,stage),...
                              Boundaries, ...
                              Mesh.Nx,Mesh.Ny);


end