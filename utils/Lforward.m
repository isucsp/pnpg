
function X=Lforward(P1,P2)
    [m2,n]=size(P1);
    [m,n1]=size(P2);
    if((n~=n1+1) || (m~=m2+1))
        error('dimensions are not consistent')
    end
    X=zeros(m,n);
    if(size(P1,1)>0)
        X(1,:)=P1(1,:); X(end,:)=-P1(end,:);
        X(2:end-1,:) = diff(P1,1,1);
    end
    Y=zeros(m,n);
    if(size(P2,2)>0)
        Y(:,1)=P2(:,1); Y(:,end)=-P2(:,end);
        Y(:,2:end-1) = diff(P2,1,2);
    end
    X=X+Y;

    % [P1;zeros(1,n)]+[P2,zeros(m,1)]-[zeros(1,n);P1]-[zeros(m,1),P2];
end

