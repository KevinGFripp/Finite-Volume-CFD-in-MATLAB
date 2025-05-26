function [uvector,vvector,P,Error] = ExplicitRungeKutta_withShape(...
                                      Solver,uvector,vvector,P, ...
                                      POperator,uOperator,vOperator, ...
                                      Mesh,Shape,Boundaries,Parameters,dt)

%% Stage 1
% % Evaluate du/dt
[Solver.du_n_dt(:,1),...
 Solver.dv_n_dt(:,1)] = NavierStokes_RK_NoP_withShape(uvector,vvector, ...
                                        uOperator,vOperator, ...
                                        Mesh,Shape,Boundaries,Parameters);
%% Stages
for S = 2:Solver.Stages
   hs = Solver.ci(S)*dt;

   uvector_s = uvector;
   vvector_s = vvector;

   for n = 1:(S-1)
       uvector_s = uvector_s + Solver.a_s(S,n)*dt*Solver.du_n_dt(:,n);
       vvector_s = vvector_s + Solver.a_s(S,n)*dt*Solver.dv_n_dt(:,n);
   end

%% Corrector
[uvector_s,vvector_s,~] = CorrectorStep_RK_withShape(uvector_s,vvector_s, ...
                                            POperator,uOperator,vOperator, ...
                                            Mesh,Shape,Boundaries,Parameters,hs);
[Solver.du_n_dt(:,S),...
 Solver.dv_n_dt(:,S)] = NavierStokes_RK_NoP_withShape(uvector_s,vvector_s, ...
                                        uOperator,vOperator, ...
                                        Mesh,Shape,Boundaries,Parameters);   
end

%% Solution
uvector = uvector + sum((dt * Solver.b.').* Solver.du_n_dt,2);
vvector = vvector + sum((dt * Solver.b.').* Solver.dv_n_dt,2);

[uvector,vvector,P] = CorrectorStep_RK_withShape(uvector,vvector, ...
                                        POperator,uOperator,vOperator, ...
                                        Mesh,Shape,Boundaries,Parameters,dt);

u_tilde = sum((dt * (Solver.b_tilde.'-Solver.b.')).* Solver.du_n_dt,2);
v_tilde = sum((dt * (Solver.b_tilde.'-Solver.b.')).* Solver.dv_n_dt,2);

u_tilde = max(abs(u_tilde));
v_tilde = max(abs(v_tilde));

Error = max(u_tilde,v_tilde);


end