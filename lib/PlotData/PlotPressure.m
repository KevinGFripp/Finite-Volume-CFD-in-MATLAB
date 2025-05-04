function PlotPressure(xscale,yscale,P_centroid,colourmap)

   imagesc(yscale,xscale,P_centroid);
   view(-90,90);
   colormap(colourmap);

end