function Lid_Driven_Cavity_withSquare()

%% -------Setup-------------
% Number of cells in grid
Nx=80;
Ny=80;

% Domain size
Width = 1.0; %m
Height =1.0; %m

% Parameters 
% Density , kinematic viscosity, time step
Parameters =struct('Nu',0.005,'rho',1,'dt',1e-4);
Reynolds_Number = 800; % Reynolds_Number = (L*u/Nu)

% Mesh Growth Rates
wx = 0.7;
wy = 0.7;

Square = DefineBody(0.5,0.2, ...
                    Width/Nx,Height/Ny, ...
                    Nx,Ny,Width/Nx,Height/Ny);

Mesh = MakeRectilinearMesh_withShape(Nx,Ny,Width,Height, ...
                                    [wx 0.75*wx wx],[wy 0.75*wy wy],Square);

% Staggered grid
uvector = zeros((Nx+1)*Ny,1); % Vx
vvector = zeros(Nx*(Ny+1),1); % Vy

% Velocity in terms of Re
velocity = Parameters.Nu * Reynolds_Number/min(Width,Height);

% Time step
max_dt = (min(Width,Height).^2)*min( (min(Mesh.u_dx)^2)/Parameters.Nu , ...
         (min(Mesh.u_dy)^2)/Parameters.Nu ); % Nu * dt / dx**2 <= 1 CFL condition
Parameters.dt = min(0.95*max_dt/log(Reynolds_Number),0.25*max_dt);

% Time to steady state
TIME = ceil(7.5*(2*Width +2*Height)/velocity/Parameters.dt);

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


 progressbar();
 for ITER = 1:TIME

 [uvector,vvector] = MakeStep_withShape(uvector,vvector,POperator,uOperator, ...
                              vOperator,Mesh,Square,Boundaries,Parameters);
  
 
   if mod(ITER,8) == 0
   u_interp_centre.Values = uvector;
   v_interp_centre.Values = vvector;

   u_centroid = reshape(u_interp_centre(Xq,Yq),Ny,Nx).';
   v_centroid = reshape(v_interp_centre(Xq,Yq),Ny,Nx).';
   figure(2)
    PlotMagnitude((1:Nx).*dx_uniform, ...
                  (1:Ny).*dy_uniform, ...
                  u_centroid,v_centroid,velocity, ...
                  colourmap);
    hold on
    PlotVectorField(Xarr,Yarr,u_centroid,v_centroid, ...
                    ceil(log10(Nx)+1),ceil(log10(Ny)+1),2.5,[0.90 0.90 1.0]);
    rectangle('Position',[Square.ycoord Square.xcoord 0.2 0.5], ...
                'FaceColor',[1 0.68 0.55],'EdgeColor','none');
    hold off
    axis off;
    set(gca,'FontSize',13,'FontName','Times','fontweight','normal');
    set(gcf,'color','w');
    drawnow()

    progressbar(ITER/TIME);
   end


 end



end


