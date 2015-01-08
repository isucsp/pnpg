function x = debias(m,n,x,A,b,M,nse)

% De-biasing for compressed sensing reconstruction
%
% Identifies the non-zero components of x and solves the reduced least-
% squares problem defined by the corresponding columns of A, assuming that 
% the number of non-zeros is in the interval (0,m].
%
% INPUTS
%
%   m, n  - A is m x n matrix; typically m < n.
%   x     - reconstructed signal x (typically sparse AND an approximate 
%           solution to Ax = b)
%   A     - m x n matrix, or function handle that implements A*x and A'*x.  
%           If the latter, this function must have the form
%
%               function y = name(x,mode)
%
%           where if mode == 1 then x is n x 1 and y = A*x, and if mode == 2 then 
%           x is m x 1 and y = A'*x.
%   b     - m x 1 vector of measurements.  Ax \approx b.
%   M     - m x m weighting matrix, unless empty.  If the former, least-
%           squares term is ||Ax - b||_M^2, if the latter, ||Ax - b||_2^2.
%   nse   - noise level for distinguishing between x's zeros and non-zeros.
%           if fpc and getM_mu are used, a multiple of sig (we use 3) is 
%           a good choice.
%
% OUTPUT
%
%   x     - components identified as approximately 0 in original signal are
%           set to exactly zero; remaining components are determined via
%           reduced least-squares problem.
%
% Copyright (C) 2007-2008 Elaine Hale, Wotao Yin and Yin Zhang
% FPC is distributed under the GNU GPL, see README.txt and COPYING.txt.

% find non-zeros and count them
inds = find(abs(x) > nse);
nnz = length(inds);

if (nnz > 0) && (nnz <= m)
    % de-bias if reduced least-squares is non-trivial and full column rank
    if ~isa(A,'function_handle')
        x = zeros(n,1);
        if isempty(M)
            x(inds) = A(:,inds)\b;
        else            
            x(inds) = (A(:,inds)'*M*A(:,inds))\(A(:,inds)'*(M*b));
        end
    else
        [xpart,flag] = lsqr(@lsqrMat,b,nse/10,20,[],[],x(inds));
        if flag < 2
            x = zeros(n,1);
            x(inds) = xpart;
        else
            warning('lsqr was not able to converge.  De-biasing has been abandoned.');
        end
    end
end

function y = lsqrMat(x,transpose)
    switch transpose
        case 'transp'
            y = A(x,2);
            y = y(inds);
        case 'notransp'
            z = zeros(n,1);
            z(inds) = x;
            y = A(z,1);
    end
end

end

% Last modified 11 January 2008.