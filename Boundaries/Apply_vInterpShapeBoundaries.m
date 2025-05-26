function [vinterp] = Apply_vInterpShapeBoundaries(vinterp,Shape,Nx,Ny)

 Left_u = index(Shape.u_grid_i(1),Shape.u_grid_j,Nx+1);
 Right_u = index(Shape.u_grid_i(end),Shape.u_grid_j,Nx+1);
 vinterp(Left_u) = 0;
 vinterp(Right_u) = 0;

end