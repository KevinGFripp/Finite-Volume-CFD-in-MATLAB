function PlotMeshandShape(Mesh,Shape)


 for n = 1:length(Mesh.xface)
     if n == 1
         hold on
     end
 plot(repmat(Mesh.xface(n),length(Mesh.yface)),Mesh.yface,'k-','LineWidth',0.25);
 end

 for m = 1:length(Mesh.yface)
 plot(Mesh.xface,repmat(Mesh.yface(m),length(Mesh.xface)),'k-','LineWidth',0.25);
 end
 
if (isempty(Shape))
else
 for y = 1:length(Shape.v_grid_j)
 plot(Mesh.xface(Shape.u_grid_i), ...
     repmat(Mesh.yface(Shape.v_grid_j(y)), ...
            length(Mesh.xface(Shape.u_grid_i))),'r-');
 end

  for x = 1:length(Shape.u_grid_i)
   plot(repmat(Mesh.xface(Shape.u_grid_i(x)), ...
               length(Shape.v_grid_j)),Mesh.yface(Shape.v_grid_j),'r-','LineWidth',0.25);
 end
 
end
 hold off

set(gca,'FontSize',14,'FontName','Times','fontweight','normal');
set(gcf,'color','w');

xlabel('x (m)');
ylabel('y (m)');

end