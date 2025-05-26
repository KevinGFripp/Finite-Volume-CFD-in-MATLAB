function [uvector,vvector,q] = CorrectorStep_RK_withShape(uvector,vvector, ...
                                      POperator,uOperator,vOperator, ...
                                      Mesh,Shape,Boundaries,Parameters,dt)

 
  [uprime,vprime] = ApplyVelocityWallBoundaries(uvector,vvector, ...
                                                Boundaries,Mesh.Nx,Mesh.Ny);
  [uprime,vprime] = ApplyVelocityShapeBoundaries(uprime,vprime, ...
                                                 Shape,Mesh.Nx,Mesh.Ny);

  % Divergence
  %% interpolate velocities to Pressure faces
  uprime_divergence = uOperator.u_interp_to_P_face * uprime;
  vprime_divergence = vOperator.v_interp_to_P_face * vprime;


  Div_uv = uOperator.Div_u_x * uprime_divergence ...
         + vOperator.Div_v_y * vprime_divergence;
  
 % Solve pressure correction
  if(Mesh.Nx*Mesh.Ny < 25000)
    q = POperator.L_pressure\(Parameters.rho*Div_uv./dt);
  else
    [q,~,~,~] = bicgstabl(POperator.L_pressure, ...
                        Parameters.rho*Div_uv./dt, ...
                        1e-8,100,POperator.L,POperator.U);
  end
 
  % update velocities
  uvector = uprime - (1/Parameters.rho) * dt * POperator.Dx_pressure * q;
  vvector = vprime - (1/Parameters.rho) * dt * POperator.Dy_pressure * q;
  
  [uvector,vvector] = ApplyVelocityWallBoundaries(uvector,vvector, ...
                                                 Boundaries,Mesh.Nx,Mesh.Ny);
  [uvector,vvector] = ApplyVelocityShapeBoundaries(uvector,vvector, ...
                                                 Shape,Mesh.Nx,Mesh.Ny);
end
