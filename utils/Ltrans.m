
function [P1,P2]=Ltrans(X)
    P1=-diff(X,1,1);
    P2=-diff(X,1,2);
    % P1=X(1:end-1,:)-X(2:end,:);
    % P2=X(:,1:end-1)-X(:,2:end);
end

