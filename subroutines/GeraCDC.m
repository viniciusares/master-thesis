function CDC=GeraCDC(nel,cdcLado1,cdcLado2,cdcLado3,cdcLado4)

CDC=zeros(4*nel,8); % col1:num, col2:tipo(code), cols3a8:cdc
for I=1:nel
    CDC(I,1)=I;
    CDC(I,2:8)=cdcLado1;
end

for I=nel+1:2*nel
    CDC(I,1)=I;
    CDC(I,2:8)=cdcLado2;
end

for I=2*nel+1:3*nel
    CDC(I,1)=I;
    CDC(I,2:8)=cdcLado3;
end

for I=3*nel+1:4*nel
    CDC(I,1)=I;
    CDC(I,2:8)=cdcLado4;
end