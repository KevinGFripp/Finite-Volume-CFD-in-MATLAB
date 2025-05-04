function Body = DefineBody(width,height,xpos,ypos,Nx,Ny,dx,dy)

xrange =((1:Nx) - (Nx/2)).*dx;
yrange =((1:Ny) - (Ny/2)).*dy;

[~,ind_i_start] = min(abs(xrange-xpos + width/2));
[~,ind_i_end] = min(abs(xrange-xpos - width/2));

[~,ind_j_start] = min(abs(yrange-ypos + height/2));
[~,ind_j_end] = min(abs(yrange-ypos - height/2));

Body = struct('grid_i',0,'grid_j',0, ...
              'u_grid_i',0,'u_grid_j',0, ...
              'v_grid_i',0,'v_grid_j',0, ...
              'xcoord',0,'ycoord',0);

Body.p_grid_i = ind_i_start:ind_i_end-1;
Body.p_grid_j = ind_j_start:ind_j_end-1;

Body.u_grid_i = ind_i_start:(ind_i_end);
% Body.u_grid_j = ind_j_start:ind_j_end;
Body.u_grid_j = Body.p_grid_j;

Body.v_grid_i = Body.p_grid_i;
% Body.v_grid_i = ind_i_start:ind_i_end;
Body.v_grid_j = ind_j_start:(ind_j_end);

Body.xcoord = xpos - width/2 +(Nx/2)*dx;
Body.ycoord = ypos - height/2 +(Ny/2)*dy;

end