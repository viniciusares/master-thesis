function Nc=Nc(ksi,tipo123)

switch tipo123
    case 1
        Nc=1/2*ksi*(ksi-1);
    case 2
        Nc=1-ksi^2;
    case 3
        Nc=1/2*ksi*(ksi+1);
end