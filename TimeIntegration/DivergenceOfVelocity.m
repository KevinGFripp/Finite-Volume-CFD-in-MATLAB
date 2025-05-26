function Div_uv = DivergenceOfVelocity(uprime,vprime,uOperator,vOperator)

  uprime_divergence = uOperator.u_interp_to_P_face * uprime;
  vprime_divergence = vOperator.v_interp_to_P_face * vprime;

  Div_uv = uOperator.Div_u_x * uprime_divergence ...
         + vOperator.Div_v_y * vprime_divergence;

end