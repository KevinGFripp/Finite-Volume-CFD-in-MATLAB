function TestInterpolation_u_Velocity_To_vgrid()
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
GrowthRatesx = 1.25;
GrowthRatesy = 1.25;

Mesh = MakeRectilinearMesh(Nx,Ny,Width,Height,GrowthRatesx,GrowthRatesy,[]);

Ru = Interpolate_u_velocity_to_v_grid(Mesh);

u_velocity = ones((Nx+1)*Ny,1).*linspace(0,1,(Nx+1)*Ny).';

figure(1)
imagesc(reshape(Ru*u_velocity,Nx,Ny+1))
view(-90,90)
set(gca,'FontSize',14,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

ExpectedSum = length(nonzeros(Ru))/4
ActualSum = full(sum(Ru,"all"))

figure(2)
imagesc(full(Ru))
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