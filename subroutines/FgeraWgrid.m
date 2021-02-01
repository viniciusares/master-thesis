function Wgrid=FgeraWgrid(w_int,xnpts,ynpts)

Wgrid=zeros(ynpts,xnpts);
for I=1:ynpts
    for J=1:xnpts
        k=xnpts*(I-1)+J;
        Wgrid(I,J)=w_int(k);
    end
end