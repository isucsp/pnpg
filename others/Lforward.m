function X=Lforward(P)

[m2,n]=size(P{1});
[m,n1]=size(P{2});

if((n~=n1+1) || (m~=m2+1))
    error('dimensions are not consistent')
end

X=[P{1};zeros(1,n)]+[P{2},zeros(m,1)]-[zeros(1,n);P{1}]-[zeros(m,1),P{2}];

end

