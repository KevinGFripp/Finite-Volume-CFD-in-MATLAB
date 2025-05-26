function [du_dt,dv_dt] = NavierStokes_RK_NoP(uvector,vvector, ...
                                             uOperator,vOperator, ...
                                             Mesh,Boundaries,Parameters)

  % u velocity 
  Del2_u = Parameters.Nu * (uOperator.L_u * uvector + uOperator.b_LU_bc);
  d_u2_dx = uvector .* (uOperator.D_ux * uvector);
  
  v_interp = Apply_vInterpWallBoundaries(vOperator.v_bar_interp * vvector, ...
                                         Boundaries,Mesh.Nx,Mesh.Ny);
 
  d_u_dy = v_interp .* (uOperator.D_uy * uvector + uOperator.b_uy_bc);

  
  
  % v velocity 
  Del2_v = Parameters.Nu * (vOperator.L_v * vvector + vOperator.b_LV_bc);
  d_v2_dy = vvector .* (vOperator.D_vy * vvector);
  
  u_interp = Apply_uInterpWallBoundaries(uOperator.u_bar_interp * uvector, ...
                                         Boundaries,Mesh.Nx,Mesh.Ny);
  
  d_v_dx = u_interp .* (vOperator.D_vx * vvector + vOperator.b_vx_bc);


  % Full Navier Stokes
  du_dt = Del2_u -d_u2_dx -d_u_dy;
  dv_dt = Del2_v -d_v2_dy -d_v_dx;

  

end
