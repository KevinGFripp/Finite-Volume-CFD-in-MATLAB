function P = InitialPressure_withShape(uvector,vvector, ...
                             POperator,uOperator,vOperator, ...
                             Mesh,Shape,Boundaries,Parameters)

[uvector,vvector] = ApplyVelocityWallBoundaries(uvector,vvector, ...
                                                Boundaries,Mesh.Nx,Mesh.Ny);
[uvector,vvector] = ApplyVelocityShapeBoundaries(uvector,vvector, ...
                                                 Shape,Mesh.Nx,Mesh.Ny);

[dudt,dvdt] = NavierStokes_RK_NoP_withShape(uvector,vvector, ...
                                    uOperator,vOperator, ...
                                    Mesh,Shape,Boundaries,Parameters);

ApplyVelocityDerivativeWallBoundaries(dudt,dvdt,Boundaries,Mesh.Nx,Mesh.Ny);

 % Divergence
  %% interpolate dudt to Pressure faces
  dudt_divergence = uOperator.u_interp_to_P_face * dudt;
  dvdt_divergence = vOperator.v_interp_to_P_face * dvdt;


  Div_uvdt = uOperator.Div_u_x * dudt_divergence ...
           + vOperator.Div_v_y * dvdt_divergence;
  
  % Solve pressure correction
  P = POperator.L_pressure\(Parameters.rho*Div_uvdt);

end