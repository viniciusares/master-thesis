% Programa apenas para surface timoshenko

xmin=-0.00;ymin=-0.00;xmax=1.0;ymax=1.0;xnpts=29;ynpts=29;

xpasso=(xmax-xmin)/(xnpts-1); ypasso=(ymax-ymin)/(ynpts-1);
xgrid=xmin:xpasso:xmax; ygrid=ymin:ypasso:ymax;
[xx,yy]=meshgrid(xgrid,ygrid);

WgridTimo=timoEALA(xx,yy,7);
surf(xx,yy,-WgridTimo)
title('esquerda:eixoY, direita:eixoX')