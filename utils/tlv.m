function out=tlv(X,type)
%This function computes the total variation of an input image X
%
% INPUT
%
% X............................. An image
% type .................... Type of total variation function. Either 'iso'
%                                  (isotropic) or 'l1' (nonisotropic)
% 
% OUTPUT
% out ....................... The total variation of X.
[m,n]=size(X);
%if(min(m,n)==1 && ~strcmpi(type,'1d'))
%    global strlen
%    fprintf('warning from tlv: using 1-d TV\n');
%    strlen=0;
%end
[P1,P2]=Ltrans(X);

switch type
    case 'iso'
        D=zeros(m,n);
        D(1:m-1,:)=P1.^2;
        D(:,1:n-1)=D(:,1:n-1)+P2.^2;
        out=sum(sum(sqrt(D)));
    case {'l1', '1d'}
        out=sum(sum(abs(P1)))+sum(sum(abs(P2)));
    otherwise
        error('Invalid total variation type. Should be either "iso" or "l1"');
end
