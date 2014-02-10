% COMPRESSED SENSING PROBLEMS
%
%--------------------------------------------------------------------------
% GENERAL DESCRIPTION AND OUTPUTS
%--------------------------------------------------------------------------
%
% function [A,b,xs,xsn,picks] = getData(m,n,k,Ameth,xmeth,varargin)
%
% Generates data for compressed sensing problems.
%
% A     - m x n matrix if A is explicit; empty otherwise
% b     - m x 1 vector equal to A*(xs+noise) + noise
% xs    - n x 1 vector with k nonzeros
% xsn   - n x 1 vector xs + noise
% picks - vector of row indices if A is implicit; empty otherwise
% 
% For compressed sensing, m <= n, and if k is small enough l1 regularized
% least squares applied to A and b will recover xs.
%
%--------------------------------------------------------------------------
% INPUTS
%--------------------------------------------------------------------------
% 
% Ameth determines how A is generated
%   0 - randn(m,n)
%   1 - randn then columns scaled to unit norm
%   2 - randn then QR to get orthonormal rows
%   3 - bernoulli +/- 1 distribution
%   4 - partial Hadamard matrix
%   5 - picks for partial fourier matrix
%   6 - picks for partial discrete cosine transform matrix 
%       REQUIRES Signal Processing Toolbox
%   7 - picks for partial 2-d discrete cosing transform matrix
%       REQUIRES Image Processing Toolbox
%
% xmeth determines how xs is generated
%   0 - randperm for support, 2*randn for values
%   1 - randperm for support, 2*(rand - 0.5) for values
%   2 - randperm for support, ones for values
%   3 - randperm for support, sign(randn) for values
%   4 - randperm for support, 1e5*randn for values
%   5 - randperm for support, 1e5*(rand - 0.5) for values
%   6 - randperm for support, 1e5 for values
%   7 - randperm for support, 1e5*sign(randn) for values
% varargin{1} = sigma1 - standard deviation of signal noise (added to xs)
% varargin{2} = sigma2 - standard deviation of meas. noise (added to b)
% varargin{3} = state - state used to initialize rand and randn (scalar
%                       integer 0 to 2^32-1)
%--------------------------------------------------------------------------

function [A,b,xs,xsn,picks] = getData(m,n,k,Ameth,xmeth,varargin)

if nargin == 8
    randn('state',varargin{3});
    rand('state',varargin{3});
end

picks = [];
switch Ameth
    case 0
        % randn, no scaling
        A = randn(m,n);
    case 1
        % randn, column scaling
        A = randn(m,n);
        for i = 1:n
            A(:,i) = A(:,i)/norm(A(:,i));
        end
    case 2
        % randn, orthonormal rows
        A = randn(m,n);
        [A,R] = qr(A',0);
        A = A';
    case 3
        % bernoulli +/- 1
        A = sign(rand([m,n]) - 0.5);
        ind = find(A == 0);
        A(ind) = ones(size(ind));
    case 4
        % partial hadamard
        A = hadamard(n);
        picks = randperm(n);
        picks = sort(picks(1:m));
        A = A(picks,:); picks = [];
    case 5
        % partial 1-d fourier transform
        picks = randperm(n);
        picks = sort(picks(1:m));
        A = A_operator(@(x) pfft(x,1,n,picks), @(x) pfft(x,2,n,picks));
    case 6
        % partial 1-d discrete cosine transform
        picks = randperm(n);
        picks = sort(picks(1:m));
        A = A_operator(@(x) pdct(x,1,n,picks), @(x) pdct(x,2,n,picks));
        %A = dctmtx(n); A = A(picks,:); 
%     case 7
%         % partial 2-d discrete cosine transform
%         picks = randperm(n);
%         picks = sort(picks(1:m));
%         A = A_operator(@(x) pdct2(x,1,n,picks), @(x) pdct2(x,2,n,picks));
end

switch Ameth
    case {0,3,4}, 
        opts_eig.tol = 1e-4; opts_eig.disp = 0; opts_eig.issym = true; opts_eig.maxit = 50;
        AAt = @(x) A*((A')*x); smax = sqrt(eigs(AAt,m,1,'lm',opts_eig));  A = (1/smax)*A;
end

% rand('state',100);
% randn('state',100);

% get original signal
p = randperm(n);
xs = zeros(n,1); kk = round(k/2);
switch xmeth
    case 0, xs(p(1:k)) = 2*randn(k,1);
    case 1, xs(p(1:k)) = 2*(rand(k,1) - 0.5);
    case 2, xs(p(1:k)) = ones(k,1);
    case 3, xs(p(1:k)) = sign(randn(k,1));
    case 4, xs(p(1:k)) = 1e5*randn(k,1);
    case 5, xs(p(1:k)) = 1e5*(rand(k,1) - 0.5);
    case 6, xs(p(1:k)) = 1e5*ones(k,1);
    case 7, xs(p(1:k)) = 1e5*sign(randn(k,1)); 
    case 8, cx = 0.5; beta = 1.5; xs = cx*[1:n]'.^(-beta);  xs(p(1:kk)) = - xs(p(1:kk)); xs(p(k+1:n)) = 0;
    case 9, cx = 1e3; beta = 1.5; xs = cx*[1:n]'.^(-beta);  xs(p(1:kk)) = - xs(p(1:kk)); xs(p(k+1:n)) = 0;
    case 10, cx = 10; beta = .5; xs = cx*exp(-[1:n]'.*beta); xs(p(1:kk)) = - xs(p(1:kk)); xs(p(k+1:n)) = 0;
    case 11, cx = 1; beta = .005; xs = cx*exp(-[1:n]'.*beta);  xs(p(1:kk)) = - xs(p(1:kk)); xs(p(k+1:n)) = 0;
    case 12, xs(p(1:kk)) = 1e5*sign(randn(kk,1)); xs(p(kk+1:2*kk)) = sign(randn(kk,1));
end


% add noise
xsn = xs;
if nargin > 5 && varargin{1} > 0
    xsn = xsn + randn(n,1)*varargin{1};
end

% get noiseless measurments
% switch Ameth
%     case {5,6,7}, b = A(false,m,n,xsn,[],picks);
%     otherwise, b = A*xsn;
% end
b = A*xsn;

% add noise
if nargin > 6 && varargin{2} > 0
    b = b + randn(m,1)*varargin{2};
end

% Omega = picks;

%save('Ameth6Xmeth2n1024m512k154seed200', 'b','Omega','n','xs');
% save('Ameth6Xmeth6n1024m512k154seed200', 'b','Omega','n','xs');

return