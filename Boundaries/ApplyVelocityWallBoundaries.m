function [uvector,vvector] = ApplyVelocityWallBoundaries(uvector,vvector, ...
                                                         Boundaries,Nx,Ny)
 Left_u = index(1,1:Ny,Nx+1);
 Right_u = index(Nx+1,1:Ny,Nx+1);

 if(strcmp(Boundaries.LeftWallType,'Wall') ...
    || strcmp(Boundaries.LeftWallType,'Inlet'))
 uvector(Left_u) = Boundaries.Left_u;
 uvector(Right_u) = Boundaries.Right_u;
 end

 Top_v = index(1:Nx,Ny+1,Nx);
 Bottom_v = index(1:Nx,1,Nx);

 vvector(Top_v) = Boundaries.Top_v;
 vvector(Bottom_v) = Boundaries.Bottom_v;

end