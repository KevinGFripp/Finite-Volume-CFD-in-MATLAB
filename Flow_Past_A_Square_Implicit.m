function Flow_Past_A_Square_Implicit()
AddLibraries();
%% -------Setup-------------
% Number of cells in grid
Nx=200;
Ny=90;

% Domain size
Width = 7.0; %m
Height =3.0; %m

% Parameters 
% Density , kinematic viscosity, time step
Parameters =struct('Nu',0.005,'rho',1,'dt',1e-4);
Reynolds_Number = 120; % Reynolds_Number = (L*u/Nu)

% Mesh Growth Rates
wx = 0.7;
wy = 0.7;

% Square Edge Length
D = 0.3;
W = 0.3;
Square = DefineBody(W,D,-2.3,Height/Ny,Nx,Ny,Width/Nx,Height/Ny);

Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height, ...
                                    [wx wx wx],[wy wy wy],Square);

Integrator = Solver('ESDIRK325',Mesh);

% Velocity in terms of Re
velocity = Parameters.Nu * Reynolds_Number/(D);

% Staggered grid
uvector = velocity*ones((Nx+1)*Ny,1); % Vx
vvector = zeros(Nx*(Ny+1),1); % Vy


% Time step
Parameters.dt = 0.5*CFL_limit(D,Mesh,Parameters);

% Stats
disp(strcat('Grid size = ',num2str(Nx),'x',num2str(Ny)));
disp(strcat(' min(dx,dy) =', ...
            num2str(EstimateMinCellSize(max(wx,wy),Nx,Width)), ...
            ' max(dx,dy) =', ...
            num2str(EstimateMaxCellSize(max(wx,wy),Nx,Width))));
disp(strcat(' max(timestep) =',num2str(CFL_limit(D,Mesh,Parameters))));
disp(strcat(' timestep =',num2str(Parameters.dt)));

% Time to steady state (s)
TIME = ceil(8.0/Parameters.dt);

% Velocity boundary values
u_left_wall = velocity;
u_right_wall = velocity;
u_top_wall = velocity;
u_bottom_wall = velocity;

v_left_wall = 0;
v_right_wall = 0;
v_top_wall = 0;
v_bottom_wall =0;

Boundaries = CreateBoundaries('Inlet',u_left_wall,v_left_wall,...
                              'Inlet',u_right_wall,v_right_wall,...
                              'Wall',u_bottom_wall,v_bottom_wall,...
                              'Wall',u_top_wall,v_top_wall);

% Matrix Operators
[POperator,uOperator,vOperator] = MakeOperators_withShape(Mesh,Boundaries,Square);
[u_xGrid,u_yGrid,Xq,Yq] = u_velocity_Interp_to_Cell_Centre_Grids(Mesh);
[v_xGrid,v_yGrid,~,~] = v_velocity_Interp_to_Cell_Centre_Grids(Mesh);

% Plotting interpolation functions
u_interp_centre = scatteredInterpolant(reshape(u_xGrid,(Nx+1)*Ny,1), ...
                                       reshape(u_yGrid,(Nx+1)*Ny,1), ...
                                       ones((Nx+1)*Ny,1));
v_interp_centre = scatteredInterpolant(reshape(v_xGrid,Nx*(Ny+1),1), ...
                                       reshape(v_yGrid,Nx*(Ny+1),1), ...
                                       ones(Nx*(Ny+1),1));
dx_uniform = Width/Nx;
dy_uniform = Height/Ny;

[Xarr,Yarr] = meshgrid((1:ceil(log10(Ny)+1):Ny).*dy_uniform, ...
                        (1:ceil(log10(Nx)+1):Nx).*dx_uniform);
colourmap = Plasma;
%% -------------------------

Pressure = InitialPressure_withShape(uvector,vvector, ...
                             POperator,uOperator,vOperator, ...
                             Mesh,Square,Boundaries,Parameters);


 times = zeros(ceil(TIME),1);
 T = 0;
 ITER = 1;
 lift =  zeros(ceil(TIME),1);
 progressbar();
 Duration = 8.0;

 while T < Duration   
 [uvector,vvector,Pressure,Error] = ImplicitRungeKutta_withShape( ...
                               Integrator,uvector,vvector,Pressure, ...
                               POperator,uOperator,vOperator, ...
                               Mesh,Square,Boundaries,Parameters,Parameters.dt);
  T = T + Parameters.dt;
  Parameters.dt = AdaptTimeStep(Error,Parameters.dt, ...
                                4*CFL_limit(D,Mesh,Parameters), ...
                                5e-3,Integrator);
  
 if(mod(ITER,1) == 0 || T == 0) 

   u_interp_centre.Values = uvector;
   v_interp_centre.Values = vvector;

   u_centroid = reshape(u_interp_centre(Xq,Yq),Ny,Nx).';
   v_centroid = reshape(v_interp_centre(Xq,Yq),Ny,Nx).';

    figure(2)
    subplot(2,1,1)
      PlotMagnitude((1:Nx).*dx_uniform,(1:Ny).*dy_uniform, ...
                  u_centroid,v_centroid,velocity,colourmap);
      clim([0 0.82*log10(Reynolds_Number)]);
      hold on
      PlotVectorField(Xarr,Yarr,u_centroid,v_centroid, ...
                     ceil(log10(Nx)+1),ceil(log10(Ny)+1),1.0,[0.1 0.3 0.1]);
%       rectangle('Position',[Square.ycoord Square.xcoord D W], ...
%                 'FaceColor',[0.3 0.99 0.3],'EdgeColor','none');
      ylim([0 6.5])
      axis off;
      hold off;
      title(strcat('Error = ',num2str(Error),' dt = ',num2str(Parameters.dt)));
      set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
      set(gcf,'color','w');

      subplot(2,1,2)
      times(ITER,1) = T;
      dat = sum(Pressure(index(Square.p_grid_i(1:end),Square.p_grid_j(1)-1,Nx)))-...
            sum(Pressure(index(Square.p_grid_i(1:end),Square.p_grid_j(end)+1,Nx)));
      lift(ITER,1) = dat;
      plot(nonzeros(times),nonzeros(lift),'k-','linewidth',1.5);
      xlim([0 Duration])
      drawnow()

    progressbar(T/Duration);
 end

ITER = ITER +1;
end


end


