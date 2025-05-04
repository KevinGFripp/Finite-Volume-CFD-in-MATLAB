function TestInterpolation_VelocityToPressureFaces()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

Width = 1.0;
Height = 1.0;
Nx = 24;
Ny = 24;
GrowthRatesx = 2;
GrowthRatesy = 2;

Mesh = MakeRectilinearMesh(Nx,Ny,Width,Height,GrowthRatesx,GrowthRatesy,[]);

%% Matlab linear interpolation
uface = interp1(Mesh.u_centre_x,...
                (Mesh.u_centre_x.*Nx).^2, ...
                Mesh.xface,"linear");
vface = interp1(Mesh.v_centre_y, ...
                (Mesh.v_centre_y.*Ny).^2, ...
                Mesh.yface,"linear");

% u velocity data
uposx = Mesh.u_centre_x;
u_velocity = zeros((Nx+1)*Ny,1);
for i = 1: Nx+1
    for j = 1:Ny
     u_velocity(index(i,j,Nx+1)) = (uposx(i).*Nx).^2;
    end
end

% v velocity data
vposy = Mesh.v_centre_y;
v_velocity = zeros(Nx*(Ny+1),1);
for i = 1: Nx
    for j = 1:Ny+1
     v_velocity(index(i,j,Nx)) = (vposy(j).*Ny).^2;
    end
end

Ru = Interpolate_u_velocity_to_face(Mesh);
Rv = Interpolate_v_velocity_to_face(Mesh);


u_interp_to_face =Ru * u_velocity;
v_interp_to_face =Rv * v_velocity;


figure(1)
plot(Mesh.xface,uface,'ko')
hold on
plot(Mesh.xface,squeeze(u_interp_to_face(index(1:Nx+1,1,Nx+1))),'r.')
hold off
set(gca,'FontSize',16,'FontName','Times','fontweight','normal');
set(gcf,'color','w');
xlabel('x face position (m)')
ylabel('(x/Width * Nx)^2')

figure(2)
plot(Mesh.yface,vface,'ko')
hold on
plot(Mesh.yface,squeeze(v_interp_to_face(index(1,1:Ny+1,Nx))),'b.')
hold off
set(gca,'FontSize',16,'FontName','Times','fontweight','normal');
set(gcf,'color','w');
xlabel('y face position (m)')
ylabel('(y/Height * Ny)^2')



end

function k = index(i,j,Nx)

k = i +(j-1)*Nx;

end