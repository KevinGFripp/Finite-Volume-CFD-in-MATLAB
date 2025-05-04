function Flow_Past_A_Square()

%% -------Setup-------------
% Number of cells in grid
Nx=180;
Ny=70;

% Domain size
Width = 6.0; %m
Height =1.5; %m

% Parameters 
% Density , kinematic viscosity, time step
Parameters =struct('Nu',0.005,'rho',1,'dt',1e-4);
Reynolds_Number = 120; % Reynolds_Number = (L*u/Nu)

% Mesh Growth Rates
wx = 0.5;
wy = 0.5;
% Square Edge Length
D = 0.25;
Square = DefineBody(D,D,-2.25,Height/Ny,Nx,Ny,Width/Nx,Height/Ny);

Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height, ...
                                    [wx wx wx],[wy wy wy],Square);

% Staggered grid
uvector = zeros((Nx+1)*Ny,1); % Vx
vvector = zeros(Nx*(Ny+1),1); % Vy

% Velocity in terms of Re
velocity = Parameters.Nu * Reynolds_Number/(D);

% Time step
max_dt = ((D).^2)*min( (min(Mesh.u_dx)^2)/Parameters.Nu , ...
         (min(Mesh.u_dy)^2)/Parameters.Nu ); % Nu * dt / dx**2 <= 1 CFL condition
Parameters.dt = min(0.95*max_dt/log(Reynolds_Number),0.25*max_dt);

% Time to steady state
TIME = ceil(7/Parameters.dt);

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

% Stats
disp(strcat('Grid size = ',num2str(Nx),'x',num2str(Ny)));
disp(strcat(' min(dx,dy) =', ...
            num2str(EstimateMinCellSize(max(wx,wy),Nx,Width)), ...
            ' max(dx,dy) =', ...
            num2str(EstimateMaxCellSize(max(wx,wy),Nx,Width))));
disp(strcat(' max(timestep) =',num2str(max_dt)));
disp(strcat(' timestep =',num2str(Parameters.dt)));

% Matrix Operators
[POperator,uOperator,vOperator] = MakeOperators_withShape(Mesh,Boundaries,Square);
[p_xGrid,p_yGrid,~,~] = P_Interp_to_Cell_Centre_Grids(Mesh);
[u_xGrid,u_yGrid,Xq,Yq] = u_velocity_Interp_to_Cell_Centre_Grids(Mesh);
[v_xGrid,v_yGrid,~,~] = v_velocity_Interp_to_Cell_Centre_Grids(Mesh);

% Plotting interpolation functions
P_interp = scatteredInterpolant(reshape(p_xGrid,(Nx)*Ny,1), ...
                                       reshape(p_yGrid,(Nx)*Ny,1), ...
                                       ones((Nx)*Ny,1));
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

  v = VideoWriter('FlowPastASquareCylinder_Re120.mp4','MPEG-4');
  v.Quality = 100;
  v.FrameRate = 60;
  open(v);

 progressbar();
 for ITER = 1:TIME

 [uvector,vvector] = MakeStep_withShape(uvector,vvector, ...
                                        POperator,uOperator,vOperator, ...
                                        Mesh,Square,Boundaries,Parameters);
  
 if(mod(ITER,7) == 0) 
   u_interp_centre.Values = uvector;
   v_interp_centre.Values = vvector;

   u_centroid = reshape(u_interp_centre(Xq,Yq),Ny,Nx).';
   v_centroid = reshape(v_interp_centre(Xq,Yq),Ny,Nx).';

    figure(2)
      PlotMagnitude((1:Nx).*dx_uniform,(1:Ny).*dy_uniform, ...
                  u_centroid,v_centroid,velocity,colourmap);
      clim([0.1 log10(Reynolds_Number)]);
      hold on
      PlotVectorField(Xarr,Yarr,u_centroid,v_centroid, ...
                     ceil(log10(Nx)+1),ceil(log10(Ny)+1),1.0,[0.1 0.1 0.1]);
      rectangle('Position',[Square.ycoord Square.xcoord 0.2 0.2], ...
                'FaceColor',[0 0.68 0.55],'EdgeColor','none');
      ylim([0 5.25])
%       axis tight;
      axis off;
      hold off;
      set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
      set(gcf,'color','w');
      set(gcf,'Position',[800 500 850 220]);
      drawnow()
      frame = getframe(gcf);
      writeVideo(v,frame);
   
    progressbar(ITER/TIME);
 end

 end
 close(v);



end


