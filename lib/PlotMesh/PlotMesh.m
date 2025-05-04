function PlotMesh(Mesh)


 for n = 1:length(Mesh.xface)
     if n == 1
         hold on
     end
 plot(repmat(Mesh.xface(n),length(Mesh.yface)),Mesh.yface,'k-','LineWidth',0.25);
 end

 for m = 1:length(Mesh.yface)
 plot(Mesh.xface,repmat(Mesh.yface(m),length(Mesh.xface)),'k-','LineWidth',0.25);
 end
 hold off

set(gca,'FontSize',14,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

xlabel('x (m)');
ylabel('y (m)');

end