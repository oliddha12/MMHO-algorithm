function drawPath(path,G)
hold on
L=size(path,1);
Sx=path(1,1)-0.5;
Sy=path(1,2)-0.5;
plot(Sx,Sy,'bs','MarkerSize',4,'LineWidth',4);   
for i=1:L-1
    figure(2)
    plot([path(i,2) path(i+1,2)]-0.5,[path(i,1) path(i+1,1)]-0.5,'Marker','+','markersize',6,'Color',[0,0.45,0.74],'LineWidth',1.5)
    hold on
end
Ex=path(end,1)-0.5;
Ey=path(end,2)-0.5;
plot(Ey,Ex,'gs','MarkerSize',4,'LineWidth',4);  
