function [du_dt,dv_dt] = NavierStokes_RK_withShape(uvector,vvector,P, ...
                                      POperator,uOperator,vOperator, ...
                                      Mesh,Shape,Boundaries,Parameters)

  % u velocity 
  Del2_u = Parameters.Nu * (uOperator.L_u * uvector + uOperator.b_LU_bc);
  d_u2_dx = uvector .* (uOperator.D_ux * uvector);
  dP_dx = (1/Parameters.rho) * POperator.Dx_pressure * P;

  v_interp = Apply_vInterpWallBoundaries(vOperator.v_bar_interp * vvector, ...
                                         Boundaries,Mesh.Nx,Mesh.Ny);
  v_interp = Apply_vInterpShapeBoundaries(v_interp,Shape,Mesh.Nx,Mesh.Ny);

  d_u_dy = v_interp .* (uOperator.D_uy * uvector + uOperator.b_uy_bc);

  
  
  % v velocity 
  Del2_v = Parameters.Nu * (vOperator.L_v * vvector + vOperator.b_LV_bc);
  d_v2_dy = vvector .* (vOperator.D_vy * vvector);
  dP_dy = (1/Parameters.rho) * POperator.Dy_pressure * P;

  u_interp = Apply_uInterpWallBoundaries(uOperator.u_bar_interp * uvector, ...
                                         Boundaries,Mesh.Nx,Mesh.Ny);
  u_interp = Apply_uInterpShapeBoundaries(u_interp,Shape,Mesh.Nx,Mesh.Ny);

  d_v_dx = u_interp .* (vOperator.D_vx * vvector + vOperator.b_vx_bc);


  % Full Navier Stokes
  du_dt = Del2_u -dP_dx -d_u2_dx -d_u_dy;
  dv_dt = Del2_v -dP_dy -d_v2_dy -d_v_dx;

  

end
