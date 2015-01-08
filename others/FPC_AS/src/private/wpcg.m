function [x,flag,relres,iter,resvec] = wpcg(A,b,tol,maxit,M1,M2,x0,varargin)

%--------------------------------------------------------------------------
% A modification of pcg in matlab
% only do one matrix-vector at each iteration
%
% Zaiwen Wen, Oct-1-2008
%--------------------------------------------------------------------------

%PCG   Preconditioned Conjugate Gradients Method.
%   X = PCG(A,B) attempts to solve the system of linear equations A*X=B for
%   X. The N-by-N coefficient matrix A must be symmetric and positive
%   definite and the right hand side column vector B must have length N.
%
%   X = PCG(AFUN,B) accepts a function handle AFUN instead of the matrix A.
%   AFUN(X) accepts a vector input X and returns the matrix-vector product
%   A*X. In all of the following syntaxes, you can replace A by AFUN.
%
%   X = PCG(A,B,TOL) specifies the tolerance of the method. If TOL is []
%   then PCG uses the default, 1e-6.
%
%   X = PCG(A,B,TOL,MAXIT) specifies the maximum number of iterations. If
%   MAXIT is [] then PCG uses the default, min(N,20).
%
%   X = PCG(A,B,TOL,MAXIT,M) and X = PCG(A,B,TOL,MAXIT,M1,M2) use symmetric
%   positive definite preconditioner M or M=M1*M2 and effectively solve the
%   system inv(M)*A*X = inv(M)*B for X. If M is [] then a preconditioner
%   is not applied. M may be a function handle MFUN returning M\X.
%
%   X = PCG(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0 is
%   [] then PCG uses the default, an all zero vector.
%
%   [X,FLAG] = PCG(A,B,...) also returns a convergence FLAG:
%    0 PCG converged to the desired tolerance TOL within MAXIT iterations
%    1 PCG iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 PCG stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during PCG became too
%      small or too large to continue computing.
%
%   [X,FLAG,RELRES] = PCG(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = PCG(A,B,...) also returns the iteration number
%   at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = PCG(A,B,...) also returns a vector of the
%   residual norms at each iteration including NORM(B-A*X0).
%
%   Example:
%      n1 = 21; A = gallery('moler',n1);  b1 = A*ones(n1,1);
%      tol = 1e-6;  maxit = 15;  M = diag([10:-1:1 1 1:10]);
%      [x1,flag1,rr1,iter1,rv1] = pcg(A,b1,tol,maxit,M);
%   Or use this parameterized matrix-vector product function:
%      afun = @(x,n)gallery('moler',n)*x;
%      n2 = 21; b2 = afun(ones(n2,1),n2);
%      [x2,flag2,rr2,iter2,rv2] = pcg(@(x)afun(x,n2),b2,tol,maxit,M);
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, QMR, SYMMLQ,
%   CHOLINC, FUNCTION_HANDLE.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.18.4.4 $ $Date: 2004/12/06 16:35:22 $

if (nargin < 2)
   error('MATLAB:pcg:NotEnoughInputs', 'Not enough input arguments.');
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
   % Check matrix and right hand side vector inputs have appropriate sizes
   [m,n] = size(A);
   if (m ~= n)
      error('MATLAB:pcg:NonSquareMatrix', 'Matrix must be square.');
   end
   if ~isequal(size(b),[m,1])
      error('MATLAB:pcg:RSHsizeMatchCoeffMatrix', ...
         ['Right hand side must be a column vector of' ...
         ' length %d to match the coefficient matrix.'],m);
   end
else
   m = size(b,1);
   n = m;
   if ~isvector(b) || (size(b,2) ~= 1) % if ~isvector(b,'column')
      error('MATLAB:pcg:RSHnotColumn',...
         'Right hand side must be a column vector.');
   end
end

% Assign default values to unspecified parameters
if (nargin < 3) || isempty(tol)
   tol = 1e-6;
end
if (nargin < 4) || isempty(maxit)
   maxit = min(n,20);
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                      % Norm of rhs vector, b
if (n2b == 0)                       % if    rhs vector is all zeros
   x = zeros(n,1);                  % then  solution is all zeros
   flag = 0;                        % a valid solution has been obtained
   relres = 0;                      % the relative residual is actually 0/0
   iter = 0;                        % no iterations need be performed
   resvec = 0;                      % resvec(1) = norm(b-A*x) = norm(0)
   if (nargout < 2)
      itermsg('pcg',tol,maxit,0,flag,iter,NaN);
   end
   return
end

if ((nargin >= 5) && ~isempty(M1))
   existM1 = 1;
   [m1type,m1fun,m1fcnstr] = iterchk(M1);
   if strcmp(m1type,'matrix')
      if ~isequal(size(M1),[m,m])
         error('MATLAB:pcg:WrongPrecondSize', ...
            ['Preconditioner must be a square matrix' ...
            ' of size %d to match the problem size.'],m);
      end
   end
else
   existM1 = 0;
   m1type = 'matrix';
end

if ((nargin >= 6) && ~isempty(M2))
   existM2 = 1;
   [m2type,m2fun,m2fcnstr] = iterchk(M2);
   if strcmp(m2type,'matrix')
      if ~isequal(size(M2),[m,m])
         error('MATLAB:pcg:WrongPrecondSize', ...
            ['Preconditioner must be a square matrix' ...
            ' of size %d to match the problem size.'],m);
      end
   end
else
   existM2 = 0;
   m2type = 'matrix';
end

if ((nargin >= 7) && ~isempty(x0))
   if ~isequal(size(x0),[n,1])
      error('MATLAB:pcg:WrongInitGuessSize', ...
         ['Initial guess must be a column vector of' ...
         ' length %d to match the problem size.'],n);
   else
      x = x0;
   end
else
   x = zeros(n,1);
end

if ((nargin > 7) && strcmp(atype,'matrix') && ...
      strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
   error('MATLAB:pcg:TooManyInputs', 'Too many input arguments.');
end

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
Ax = iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
r = b - Ax;
normr = norm(r);                   % Norm of residual

if (normr <= tolb)                 % Initial guess is a good enough solution
   flag = 0;
   relres = normr / n2b;
   iter = 0;
   resvec = normr;
   if (nargout < 2)
      itermsg('pcg',tol,maxit,0,flag,iter,relres);
   end
   return
end

resvec = zeros(maxit+1,1);         % Preallocate vector for norm of residuals
resvec(1) = normr;                 % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of minimum residual
rho = 1;
stag = 0;                          % stagnation of the method

% loop over maxit iterations (unless convergence or failure)

for i = 1 : maxit
   Axp = Ax; % save A*x
   if existM1
      y = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
      if any(~isfinite(y))
         flag = 2;
         break
      end
   else % no preconditioner
      y = r;
   end

   if existM2
      z = iterapp('mldivide',m2fun,m2type,m2fcnstr,y,varargin{:});
      if any(~isfinite(z))
         flag = 2;
         break
      end
   else % no preconditioner
      z = y;
   end

   rho1 = rho;
   rho = r' * z;
   if ((rho == 0) || isinf(rho))
      flag = 4;
      break
   end
   if (i == 1)
      p = z;
   else
      beta = rho / rho1;
      if ((beta == 0) || isinf(beta))
         flag = 4;
         break
      end
      p = z + beta * p;
   end
   q = iterapp('mtimes',afun,atype,afcnstr,p,varargin{:}); % q = A*p
   pq = p' * q;
   if ((pq <= 0) || isinf(pq))
      flag = 4;
      break
   else
      alpha = rho / pq;
   end
   if isinf(alpha)
      flag = 4;
      break
   end
   if (alpha == 0)                  % stagnation of the method
      stag = 1;
   end

   % Check for stagnation of the method
   if (stag == 0)
      stagtest = zeros(n,1);
      ind = (x ~= 0);
      stagtest(ind) = p(ind) ./ x(ind);
      stagtest(~ind & p ~= 0) = Inf;
      if (abs(alpha)*norm(stagtest,inf) < eps)
         stag = 1;
      end
   end

   x = x + alpha * p;               % form new iterate
   Ax = Axp + alpha*q;
   %normr = norm(b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:}));
   normr = norm(b - Ax);
   resvec(i+1) = normr;

   if (normr <= tolb)               % check for convergence
      flag = 0;
      iter = i;
      break
   end

   if (stag == 1)
      flag = 3;
      break
   end

   if (normr < normrmin)           % update minimal norm quantities
      normrmin = normr;
      xmin = x;
      imin = i;
   end

   r = r - alpha * q;

end                                % for i = 1 : maxit

% returned solution is first with minimal residual
if (flag == 0)
   relres = normr / n2b;
else
   x = xmin;
   iter = imin;
   relres = normrmin / n2b;
end

% truncate the zeros from resvec
if ((flag <= 1) || (flag == 3))
   resvec = resvec(1:i+1);
else
   resvec = resvec(1:i);
end

% only display a message if the output flag is not used
if (nargout < 2)
   itermsg('pcg',tol,maxit,i,flag,iter,relres);
end

end

function [atype,afun,afcnstr] = iterchk(A)
%ITERCHK  Checks arguments to iterative methods.
%   [ATYPE,AFUN,AFCNSTR] = ITERCHK(A) returns the following:
%   ATYPE is either 'matrix', 'function', 'expression' or 'inline object'.
%   AFUN is the function name or inline object.
%   AFUN is '' if ATYPE is 'matrix'.
%   AFCNSTR is the function name if ATYPE is 'function'.
%   AFCNSTR is the formula of the function if ATYPE is 'expression' or
%   'inline object'.  AFCNSTR is '' if ATYPE is 'matrix'.
%
%   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 1.8.4.2 $ $Date: 2004/12/06 16:35:56 $


[afun,afunmsg] = fcnchk(A);
if isempty(afunmsg)
   if isa(afun,'inline')      
      if isa(A,'inline')
         atype = 'inline object';
      else
         atype = 'expression';
      end
      afcnstr = formula(afun);
   else % both function_handles @fun and function names 'fun'
      atype = 'function';
      if isa(A,'function_handle')
          afcnstr = func2str(A);
      else
          afcnstr = A;
      end
   end
elseif isa(A,'float')
   afun = A;
   atype = 'matrix';
   afcnstr = '';
else
   error('MATLAB:iterchk:InvalidInput',...
         'Argument must be a floating point matrix or a function handle.');
end

end


function y = iterapp(op,afun,atype,afcnstr,x,varargin)
%ITERAPP   Apply matrix operator to vector and error gracefully.
%   ITERAPP(OP,AFUN,ATYPE,AFCNSTR,X) applies matrix operator AFUN to vector
%   X. If ATYPE is 'matrix, then AFUN is a matrix and the OP is applied
%   directly. OP is either 'mtimes' or 'mldivide'.
%   ATYPE and AFCNSTR are used in case of error.
%   ITERAPP(OP,AFUN,ATYPE,AFCNSTR,X,P1,P2,...) allows extra arguments to
%   AFUN(X,P1,P2,...) although this usage is now discouraged in favor of
%   using anonymous functions.
%   AFUN(X,P1,P2,...,PN,TFLAG) should accept a TFLAG as its final input if
%   the calling function is BICG, LSQR or QMR. TFLAG is either 'transp' or
%   'notransp' depending on whether A' OP X or A OP X is required.
%   ITERAPP is designed for use by iterative methods like PCG which
%   require matrix operators AFUN representing matrices A to operate on
%   vectors X and return A*X and may also have operators MFUN representing
%   preconditioning matrices M operate on vectors X and return M\X.
%
%   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ.

%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.7.4.3 $ $Date: 2005/06/21 19:37:27 $

if strcmp(atype,'matrix')
    switch lower(op)
        case 'mtimes'
            if (nargin >= 6) && isequal(varargin{end},'transp')
                y = afun' * x;
            else
                y = afun * x;
            end
        case 'mldivide'
            if (nargin >= 6) && isequal(varargin{end},'transp')
                y = afun' \ x;
            else
                y = afun \ x;
            end
        otherwise
            error('MATLAB:iterapp:InvalidOp', 'Invalid operation.')
    end
else
    try
        if (nargin >= 6) && isequal(varargin{end},'notransp')
            % A request for A*x coming from BICG, LSQR and QMR
            try
                % New syntax: we now request afun(x,P1,P2,...,PN,'notransp')
                y = afun(x,varargin{:});
            catch
                % Old syntax: we used to accept afun(x,P1,P2,...,PN)
                y = afun(x,varargin{1:end-1});
            end
        else
            % A request for A*x
            % coming from BICGSTAB, CGS, GMRES, MINRES, PCG or SYMMLQ
            % with the call always afun(P1,P2,...,PN)
            % or a request for A'*x coming from
            % BICG, LSQR and QMR in the afun(x,P1,P2,...,PN,'transp') case
            y = afun(x,varargin{:});
        end
    catch
        error('MATLAB:InvalidInput', ['user supplied %s ==> %s\n' ...
            'failed with the following error:\n\n%s'], ...
            atype,afcnstr,lasterr);
    end

    if ~isvector(y) || (size(y,2) ~= 1)
        error('MATLAB:MustReturnColumn', ['user supplied %s ==> %s\n' ...
            'must return a column vector.'], ...
            atype,afcnstr)
    end
end

end
