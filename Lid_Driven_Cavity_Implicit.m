function Lid_Driven_Cavity_Implicit()
AddLibraries();
%% -------Setup-------------
% Number of cells in grid
Nx=70;
Ny=70;

% Domain size
Width = 1.0; %m
Height =1.0; %m

% Parameters 
% Density , kinematic viscosity, time step
Parameters =struct('Nu',0.005,'rho',1,'dt',1e-4);
Reynolds_Number = 1000; % Reynolds_Number = (L*u/Nu)

% Mesh Growth Rates
wx = 0.9;
wy = 0.9;
Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height,wx,wy,[]);

Integrator = Solver('ESDIRK325',Mesh);

% Staggered grid
uvector = zeros((Nx+1)*Ny,1); % Vx
vvector = zeros(Nx*(Ny+1),1); % Vy

% Velocity in terms of Re
velocity = Parameters.Nu * Reynolds_Number/min(Width,Height);

% Time step
Parameters.dt = 0.5*Fixed_dt(CFL_limit(Width,Mesh,Parameters),Reynolds_Number);

% Stats
disp(strcat('Grid size = ',num2str(Nx),'x',num2str(Ny)));
disp(strcat(' min(dx,dy) =', ...
            num2str(EstimateMinCellSize(max(wx,wy),Nx,Width)), ...
            ' max(dx,dy) =', ...
            num2str(EstimateMaxCellSize(max(wx,wy),Nx,Width))));
disp(strcat(' max(timestep) =',num2str(CFL_limit(Width,Mesh,Parameters))));
disp(strcat(' timestep =',num2str(Parameters.dt)));

% Velocity boundary values
u_left_wall = 0;
u_right_wall = 0;
u_top_wall = velocity;
u_bottom_wall = 0;

v_left_wall = 0;
v_right_wall = 0;
v_top_wall = 0;
v_bottom_wall =0;

Boundaries = CreateBoundaries('Wall',u_left_wall,v_left_wall,...
                              'Wall',u_right_wall,v_right_wall,...
                              'Wall',u_bottom_wall,v_bottom_wall,...
                              'Inlet',u_top_wall,v_top_wall);

% Matrix Operators
[POperator,uOperator,vOperator] = MakeOperators(Mesh,Boundaries);
[u_xGrid,u_yGrid,Xq,Yq] = u_velocity_Interp_to_Cell_Centre_Grids(Mesh);
[v_xGrid,v_yGrid,~,~] = v_velocity_Interp_to_Cell_Centre_Grids(Mesh);

% Plotting interpolation functions
u_interp_centre = scatteredInterpolant(u_xGrid(:),u_yGrid(:), ...
                                        ones((Nx+1)*Ny,1));
v_interp_centre = scatteredInterpolant(v_xGrid(:),v_yGrid(:), ...
                                        ones(Nx*(Ny+1),1));
dx_uniform = Width/Nx;
dy_uniform = Height/Ny;

[Xarr,Yarr] = meshgrid((1:ceil(log10(Ny)+1):Ny).*dy_uniform, ...
                        (1:ceil(log10(Nx)+1):Nx).*dx_uniform);
colourmap = Plasma;
%% -------------------------

 D1 = load('Ghia_LidDrivenCavity_u_velocity_x=0.5.txt');
 D2 = load('Ghia_LidDrivenCavity_v_velocity_y=0.5.txt');
 Analytic_y = D1(:,1);
 Analytic_x = D2(:,1);
 Analytic_u = D1(:,4);
 Analytic_v = D2(:,4);
 [~,Indx] = min(abs(Mesh.u_centre_x - 0.5*Width));
 [~,Indy] = min(abs(Mesh.v_centre_y - 0.5*Height));

 Pressure = InitialPressure(uvector,vvector, ...
                            POperator,uOperator,vOperator, ...
                            Mesh,Boundaries,Parameters);

 ITER = 1;
 T = 0;
 Duration = (8.0*(2*Width +2*Height)/velocity);
 progressbar();
 while T < Duration
 
 [uvector,vvector,Pressure,Error] =ImplicitRungeKutta_NoShape( ...
                               Integrator,uvector,vvector,Pressure, ...
                               POperator,uOperator,vOperator, ...
                               Mesh,Boundaries,Parameters,Parameters.dt);

   T = T + Parameters.dt;

   Parameters.dt = AdaptTimeStep(Error,Parameters.dt, ...
                                 2*CFL_limit(Width,Mesh,Parameters), ...
                                 1e-2,Integrator);

  if (mod(ITER,1) == 0 || T == 0)
   
   % interp velocity to Pressure centres
   u_interp_centre.Values = uvector;
   v_interp_centre.Values = vvector;

   u_centroid = reshape(u_interp_centre(Xq,Yq),Ny,Nx).';
   v_centroid = reshape(v_interp_centre(Xq,Yq),Ny,Nx).';
 
   figure(2)
   subplot(1,2,1)
    PlotMagnitude((1:Nx).*dx_uniform, ...
                  (1:Ny).*dy_uniform, ...
                  u_centroid,v_centroid,velocity, ...
                  colourmap);
    hold on
    PlotVectorField(Xarr,Yarr,u_centroid,v_centroid, ...
                    ceil(log10(Nx)+1),ceil(log10(Ny)+1),2.5,[0.90 0.90 1.0]);
    hold off
    axis off;
    set(gcf,'color','w');
    title(strcat('Error = ',num2str(Error),' dt = ',num2str(Parameters.dt)));

    subplot(1,2,2)
     plot(Mesh.u_centre_y, ...
          squeeze(uvector(index(Indx,1:Ny,Nx+1)))./velocity, ...
          'b-','linewidth',2);
     hold on
     plot(Mesh.v_centre_x, ...
          squeeze(vvector(index(1:Nx,Indy,Nx)))./velocity, ...
          'r-','linewidth',2);
     plot(Analytic_x,Analytic_v,'r','Marker','^', ...
       'MarkerSize',8,'linestyle','none','MarkerFaceColor','r');
     plot(Analytic_y,Analytic_u,'b','Marker','>', ...
       'MarkerSize',8,'linestyle','none','MarkerFaceColor','b');
     hold off
     xlabel('y/Height')
     ylabel('u/max(u)')
     set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
     set(gcf,'color','w');
     drawnow();

   progressbar(T/Duration);

  end
  ITER = ITER +1;
 end

end


