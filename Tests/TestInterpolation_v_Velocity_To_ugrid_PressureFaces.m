function TestInterpolation_v_Velocity_To_ugrid_PressureFaces()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

Width = 1.0;
Height = 1.0;
Nx = 16;
Ny = 16;
GrowthRatesx = 2;
GrowthRatesy = 2;

Mesh = MakeRectilinearMesh(Nx,Ny,Width,Height,GrowthRatesx,GrowthRatesy,[]);

Rv = Interpolate_v_velocity_to_face(Mesh);

v_velocity = ones(Nx*(Ny+1),1).*linspace(0,1,Nx*(Ny+1)).';

figure(1)
imagesc(reshape(Rv*v_velocity,Nx+1,Ny))
view(-90,90)
set(gca,'FontSize',14,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

sum(Rv,"all")
length(nonzeros(Rv))
figure(2)
imagesc(full(Rv))
set(gca,'FontSize',14,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

ExpectedSum = length(nonzeros(Rv))/2 +3/2 * Ny
ActualSum = full(sum(Rv,"all"))


end

function k = index(i,j,Nx)

k = i +(j-1)*Nx;

end