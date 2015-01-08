function [mu,A,b] = scaleProblem(n,mu,A,b,M)

% Scales mu, A and b so that the largest eigenvalue of A'*M*A is 1 and the
% new problem
%
%       min ||x||_1 + mu/2 * ||Ax - b||_M^2
%
% is equivalent to the old one.  We assume that A and M are both functions
% that take the arguments (x,mode), and return, for instance, A*x if mode
% == 1 and A'*x if mode == 2.

eopts.disp = 0;
eopts.issym = true;
if ~isreal(A(rand(n,1),1))
    eopts.isreal = false;
end

if isempty(M)
    fh = @(x) A(A(x,1),2);
else
    fh = @(x) A(M(A(x,1),1),2);
end
s2 = eigs(fh,n,1,'lm',eopts);
if s2 > 1 - eps
    mu = mu*s2;
    b = b/sqrt(s2);
    A = @(x,mode) A(x,mode)/sqrt(s2);
end

return