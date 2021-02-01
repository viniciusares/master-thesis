function JACxy=JACxy(X1,X2,X3,ksi,xy)

x1=X1(1);y1=X1(2);
x2=X2(1);y2=X2(2);
x3=X3(1);y3=X3(2);

dN1c=(2*ksi-1)/2;
dN2c=-2*ksi;
dN3c=(2*ksi+1)/2;
dxdksi=x1*dN1c+x2*dN2c+x3*dN3c;
dydksi=y1*dN1c+y2*dN2c+y3*dN3c;
switch xy
    case 0
        JACxy=sqrt(dxdksi^2+dydksi^2);
    case 1
        JACxy=dxdksi;
    case 2
        JACxy=dydksi;
end