function PlotMagnitude(xscale,yscale,u_centroid,v_centroid,velocity,colourmap)

   Energy = u_centroid.^2 + v_centroid.^2;
   imagesc(yscale,xscale,sqrt(Energy)./velocity);
   view(-90,90);
   colormap(colourmap);
   clim([0 1]);

end