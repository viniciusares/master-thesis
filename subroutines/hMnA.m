function hMnA=hMnA(tipo246,X1,X3,v,k0)

x1=X1(1); y1=X1(2); x3=X3(1); y3=X3(2);
L=sqrt((x3-x1)^2+(y3-y1)^2);
v1=v+1;

switch tipo246
    case 2
        %h12
        hMnA=L/(64*pi)*(...
        -6*v1*k0^2-3*v1*(k0-1)*k0^2*log(1-k0)+...
            3*v1*((k0-1)*k0^2+2)*log(k0+1)+6*v1*k0+7*v+6*v1*log(L/2)+1);
    case 4
        %h14
        hMnA=L/(32*pi)*(...
            6*v1*k0^2+v*log(1/4*(1-k0^2))-2*v1*(3*k0^2-4)*k0*atanh(k0)-...
            3*v+log(-1/4*(k0-1)*(k0+1))+2*v1*log(L)-5);
    case 6      
        %h16
        hMnA=L/(64*pi)*(...
            6*v1*(k0+1)*k0^2*atanh(k0)+...
            v*(7-6*k0*(k0+1))-6*(k0+1)*k0+6*v1*log(1/2*(L-k0*L))+1);
end