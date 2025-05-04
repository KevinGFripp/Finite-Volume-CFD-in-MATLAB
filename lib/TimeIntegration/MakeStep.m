function [uvector,vvector] = MakeStep(uvector,vvector, ...
                                      POperator,uOperator,vOperator, ...
                                      Mesh,Boundaries,Parameters)

  % u velocity 
  Del2_u = Parameters.Nu * (uOperator.L_u * uvector + uOperator.b_LU_bc);
  d_u2_dx = uvector .* (uOperator.D_ux * uvector);
  v_interp = Apply_vInterpWallBoundaries(vOperator.v_bar_interp * vvector, ...
                                         Boundaries,Mesh.Nx,Mesh.Ny);
  d_u_dy = v_interp .* (uOperator.D_uy * uvector + uOperator.b_uy_bc);

  uprime = uvector + Parameters.dt * (Del2_u -d_u2_dx -d_u_dy);


  % v velocity 
  Del2_v = Parameters.Nu * (vOperator.L_v * vvector + vOperator.b_LV_bc);
  d_v2_dy = vvector .* (vOperator.D_vy * vvector);
  u_interp = Apply_uInterpWallBoundaries(uOperator.u_bar_interp * uvector, ...
                                         Boundaries,Mesh.Nx,Mesh.Ny);
  d_v_dx = u_interp .* (vOperator.D_vx * vvector + vOperator.b_vx_bc);

  vprime = vvector + Parameters.dt * (Del2_v -d_v2_dy -d_v_dx);

  [uprime,vprime] = ApplyVelocityWallBoundaries(uprime,vprime, ...
                                                Boundaries,Mesh.Nx,Mesh.Ny);

  % Divergence
  %% interpolate velocities to Pressure faces
  uprime_divergence = uOperator.u_interp_to_P_face * uprime;
  vprime_divergence = vOperator.v_interp_to_P_face * vprime;

  Div_uv = uOperator.Div_u_x * uprime_divergence ...
         + vOperator.Div_v_y * vprime_divergence;
  
  % Solve pressure correction
  q = POperator.L_pressure\(Parameters.rho*Div_uv./Parameters.dt);

  % update velocities
  uvector = uprime - (1/Parameters.rho) * Parameters.dt * POperator.Dx_pressure * q;
  vvector = vprime - (1/Parameters.rho) * Parameters.dt * POperator.Dy_pressure * q;
  [uvector,vvector] = ApplyVelocityWallBoundaries(uvector,vvector, ...
                                                 Boundaries,Mesh.Nx,Mesh.Ny);


end