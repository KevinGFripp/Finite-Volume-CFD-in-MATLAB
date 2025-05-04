function [vinterp] = Apply_vInterpWallBoundaries(vinterp,Boundaries,Nx,Ny)

 Left_u = index(1,1:Ny,Nx+1);
 Right_u = index(Nx+1,1:Ny,Nx+1);

 if(strcmp(Boundaries.LeftWallType,'Wall') ...
    || strcmp(Boundaries.LeftWallType,'Inlet'))
 vinterp(Left_u) = Boundaries.Left_v;
 end

 if(strcmp(Boundaries.RightWallType,'Wall') ...
    || strcmp(Boundaries.RightWallType,'Inlet'))
 vinterp(Right_u) = Boundaries.Right_v;
 end

end