function CDC=GeraCDCxy(nelX,nelY,cdcLado1,cdcLado2,cdcLado3,cdcLado4)

CDC=zeros(2*nelX+2*nelY,8); % col1:num, col2:tipo(code), cols3a8:cdc
for I=1:nelX
    CDC(I,1)=I;
    CDC(I,2:8)=cdcLado1;
end

for I=nelX+1:nelX+nelY
    CDC(I,1)=I;
    CDC(I,2:8)=cdcLado2;
end

for I=nelX+nelY+1:nelX+nelY+nelX
    CDC(I,1)=I;
    CDC(I,2:8)=cdcLado3;
end

for I=nelX+nelY+nelX+1:nelX+nelY+nelX+nelY
    CDC(I,1)=I;
    CDC(I,2:8)=cdcLado4;
end