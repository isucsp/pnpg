function P=Ltrans(X)
P{1}=X(1:end-1,:)-X(2:end,:);
P{2}=X(:,1:end-1)-X(:,2:end);
