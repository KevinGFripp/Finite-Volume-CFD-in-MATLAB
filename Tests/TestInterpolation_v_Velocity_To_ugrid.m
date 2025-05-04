function TestInterpolation_v_Velocity_To_ugrid()
%% include lib
basefolder = ...
['C:\Users\kevin\Documents\MANNGA\Data\April 25\' ...
'2D_StaggeredGrid_NonUniformMesh_FiniteVolume\'];
addpath(genpath(basefolder));
%%

Width = 1.0;
Height = 1.0;
Nx = 32;
Ny = 32;
GrowthRatesx = 1.25;
GrowthRatesy = 1.25;

Mesh = MakeRectilinearMesh(Nx,Ny,Width,Height,GrowthRatesx,GrowthRatesy,[]);

Rv = Interpolate_v_velocity_to_u_grid(Mesh);

v_velocity = ones(Nx*(Ny+1),1).*linspace(0,1,Nx*(Ny+1)).';

figure(1)
imagesc(reshape(Rv*v_velocity,Nx+1,Ny))
view(-90,90)
set(gca,'FontSize',14,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

ExpectedSum = length(nonzeros(Rv))/4
ActualSum = full(sum(Rv,"all"))

figure(2)
imagesc(full(Rv))
set(gca,'FontSize',14,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

% figure(3)
% PlotMesh(Mesh)
% hold on
% PlotMeshCentres(Mesh)
% hold off



end

function k = index(i,j,Nx)

k = i +(j-1)*Nx;

end