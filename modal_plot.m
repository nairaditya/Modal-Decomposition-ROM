function modal_plot(geom, x_grid,y_grid,Mode)
%% plotting routine for modes
figure;
grey = [0.9,0.9,0.9];
hh1 = pcolor(x_grid,y_grid,real(Mode));
colormap(jet);
set(hh1,'Edgecolor','none');
hold on;
fill(geom(:,1),geom(:,2),grey);
xlim([-0.3 2]);
ylim([-0.5 0.5]);
axis equal;
axis off;
set(gca,'position',[0 0 1 1],'units','normalized')
