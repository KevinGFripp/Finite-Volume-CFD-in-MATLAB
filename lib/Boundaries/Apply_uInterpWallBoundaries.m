function [uinterp] = Apply_uInterpWallBoundaries(uinterp,Boundaries,Nx,Ny)

 TopWall_v = index(1:Nx,Ny+1,Nx);
 BottomWall_v = index(1:Nx,1,Nx);
 uinterp(TopWall_v) = Boundaries.Top_u;
 uinterp(BottomWall_v) = Boundaries.Bottom_u;

end