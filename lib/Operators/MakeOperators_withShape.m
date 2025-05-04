function [POperator,uOperator,vOperator] = ...
         MakeOperators_withShape(Mesh,Boundaries,Shape)

%% Matrix Operators

% Pressure
L_pressure = PressureLaplacian_withShape(Mesh,Shape);
Dx_pressure = Pressure_Gradientx_CellWise_withShape(Mesh,Shape);
Dy_pressure = Pressure_Gradienty_CellWise_withShape(Mesh,Shape);

POperator.L_pressure = L_pressure;
POperator.Dx_pressure = Dx_pressure;
POperator.Dy_pressure = Dy_pressure;

% u velocity 
[L_u,b_LU_bc] = uVelocityLaplacian_withShape(Mesh,Boundaries,Shape);
D_ux = uVelocityGradient_x_withShape(Mesh,Shape);
[D_uy,b_uy_bc] = uVelocityGradient_y_withShape(Mesh,Boundaries,Shape);
u_bar_interp = Interpolate_u_velocity_to_v_grid(Mesh);
Div_u_x = Divergence_x_withShape(Mesh,Shape);
u_interp_to_P_face = Interp_u_velocity_to_face_withShape(Mesh,Shape);

uOperator.L_u = L_u;
uOperator.b_LU_bc = b_LU_bc;
uOperator.D_ux = D_ux;
uOperator.D_uy = D_uy;
uOperator.b_uy_bc = b_uy_bc;
uOperator.u_bar_interp = u_bar_interp;
uOperator.Div_u_x =Div_u_x;
uOperator.u_interp_to_P_face = u_interp_to_P_face;


% v velocity
[L_v,b_LV_bc] = vVelocityLaplacian_withShape(Mesh,Boundaries,Shape);
[D_vx,b_vx_bc] = vVelocityGradient_x_withShape(Mesh,Boundaries,Shape);
D_vy = vVelocityGradient_y_withShape(Mesh,Shape);
v_bar_interp = Interpolate_v_velocity_to_u_grid(Mesh);
Div_v_y = Divergence_y_withShape(Mesh,Shape);
v_interp_to_P_face = Interp_v_velocity_to_face_withShape(Mesh,Shape);

vOperator.L_v = L_v;
vOperator.b_LV_bc = b_LV_bc;
vOperator.D_vx = D_vx;
vOperator.b_vx_bc = b_vx_bc;
vOperator.D_vy = D_vy;
vOperator.v_bar_interp = v_bar_interp;
vOperator.Div_v_y =Div_v_y;
vOperator.v_interp_to_P_face = v_interp_to_P_face;


%% -------------------------


end