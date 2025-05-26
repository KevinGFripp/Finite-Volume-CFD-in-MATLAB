function [uinterp] = Apply_uInterpShapeBoundaries(uinterp,Shape,Nx,Ny)

 TopWall_v = index(Shape.v_grid_i,Shape.v_grid_j(end),Nx);
 BottomWall_v = index(Shape.v_grid_i,Shape.v_grid_j(1),Nx);
 uinterp(TopWall_v) = 0;
 uinterp(BottomWall_v) = 0;

end