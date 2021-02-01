function Nd=Nd(ksi,tipo123)

switch tipo123
    case 1
        Nd=ksi*(ksi*9/8-3/4);
    case 2
        Nd=(1-ksi*3/2)*(1+ksi*3/2);
    case 3
        Nd=ksi*(ksi*9/8+3/4);
end