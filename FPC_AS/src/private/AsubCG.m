function Asubx = AsubCG(x)    
global MT AT nfe_sub nge_sub;

if isstruct(AT);
    xtmp=zeros(AT.n,1); xtmp(AT.col)=x; Asubx1 = AT.A*xtmp;
    if isempty(MT); Asubx = (AT.A')*Asubx1;  else Asubx = (AT.A')*(MT*Asubx1);  end;    Asubx = Asubx(AT.col);
else
    Asubx1 = AT*x; if isempty(MT);  Asubx = AT'*Asubx1;   else Asubx = AT'*(MT*Asubx1); end
end
nfe_sub = nfe_sub + 1;  nge_sub = nge_sub + 1;
end