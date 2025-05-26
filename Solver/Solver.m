function S = Solver(type,Mesh)

S = struct('Name',type,'Stages',1,'du_n_dt',0,'dv_n_dt',0,'Order',1);

if strcmp(type,'RKBS23')
   S.Stages = 4;
   S.Order = 3;
   [S.ci,S.a_s,S.b,S.b_tilde] = BogackiShampine32_Coefficients();
   S.LastStageSolution = true;
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);
end

if strcmp(type,'RK45Fehlberg')
   S.Stages = 6;
   S.Order = 5;
   [S.ci,S.a_s,S.b,S.b_tilde] = RK45Fehlberg_Coefficients();
   S.LastStageSolution = false;
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);
end

if strcmp(type,'ESDIRK325')
   S.Stages = 5;
   S.Order = 3;
   [S.ci,S.a_s,S.b,S.b_tilde] = ESDIRK325_Coefficients();
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);

end

if strcmp(type,'ESDIRK32')
   S.Stages = 4;
   S.Order = 3;
   [S.ci,S.a_s,S.b,S.b_tilde] = ESDIRK32_Coefficients();
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);

end

if strcmp(type,'ESDIRK32c')
   S.Stages = 4;
   S.Order = 3;
   [S.ci,S.a_s,S.b,S.b_tilde] = ESDIRK32c_Coefficients();
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);

end

if strcmp(type,'ESDIRK54a')
   S.Stages = 7;
   S.Order = 5;
   [S.ci,S.a_s,S.b,S.b_tilde] = ESDIRK54a_Coefficients();
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);
end

if strcmp(type,'ESDIRK547L')
   S.Stages = 7;
   S.Order = 5;
   [S.ci,S.a_s,S.b,S.b_tilde] = ESDIRK547L_Coefficients();
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);
end

if strcmp(type,'ESDIRK548L')
   S.Stages = 8;
   S.Order = 5;
   [S.ci,S.a_s,S.b,S.b_tilde] = ESDIRK548L_Coefficients();
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);
end

if strcmp(type,'ESDIRK436L')
   S.Stages = 6;
   S.Order = 4;
   [S.ci,S.a_s,S.b,S.b_tilde] = ESDIRK436L_Coefficients();
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);
end

if strcmp(type,'ESDIRK437L')
   S.Stages = 7;
   S.Order = 4;
   [S.ci,S.a_s,S.b,S.b_tilde] = ESDIRK437L_Coefficients();
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);
end

if strcmp(type,'RK4')
   S.Stages = 4;
   S.Order = 4;
   S.LastStageSolution = false;
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);
end

if strcmp(type,'Euler')
   S.Stages = 1;
   S.Order = 1;
   S.LastStageSolution = true;
   S.du_n_dt = zeros((Mesh.Nx+1)*Mesh.Ny,S.Stages);
   S.dv_n_dt = zeros(Mesh.Nx*(Mesh.Ny+1),S.Stages);
end

end