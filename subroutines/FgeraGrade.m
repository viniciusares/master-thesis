function P_int=FgeraGrade(xmin,ymin,xmax,ymax,xnpts,ynpts)

cont=0; P_int=zeros(xnpts*ynpts,3);
for I=1:ynpts % loop in y
    for J=1:xnpts % loop in x
        cont=cont+1;
        x=xmin+(J-1)*(xmax-xmin)/(xnpts-1);
        y=ymin+(I-1)*(ymax-ymin)/(ynpts-1);
        P_int(cont,:)=[cont,x,y];
    end
end