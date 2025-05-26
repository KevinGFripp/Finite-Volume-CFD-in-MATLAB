function [dudt,dvdt] = ApplyVelocityDerivativeShapeBoundaries(dudt,dvdt, ...
                                                          Shape,Nx,Ny)

 Left_u = index(Shape.u_grid_i(1),Shape.u_grid_j,Nx+1);
 Right_u = index(Shape.u_grid_i(end),Shape.u_grid_j,Nx+1);

 dudt(Left_u) = 0;
 dudt(Right_u) = 0;

 Top_v = index(Shape.v_grid_i,Shape.v_grid_j(end),Nx);
 Bottom_v = index(Shape.v_grid_i,Shape.v_grid_j(1),Nx);

 dvdt(Top_v) = 0;
 dvdt(Bottom_v) = 0;

end