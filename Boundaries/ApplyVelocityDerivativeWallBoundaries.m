function [dudt,dvdt] = ApplyVelocityDerivativeWallBoundaries( ...
                                                         dudt,dvdt, ...
                                                         Boundaries,Nx,Ny)
 Left_u = index(1,1:Ny,Nx+1);
 Right_u = index(Nx+1,1:Ny,Nx+1);

 if(strcmp(Boundaries.LeftWallType,'Wall') ...
    || strcmp(Boundaries.LeftWallType,'Inlet'))
 dudt(Left_u) = 0;
 dudt(Right_u) = 0;
 end

 Top_v = index(1:Nx,Ny+1,Nx);
 Bottom_v = index(1:Nx,1,Nx);

 dvdt(Top_v) = 0;
 dvdt(Bottom_v) = 0;

end