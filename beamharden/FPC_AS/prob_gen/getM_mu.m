% RECOMMENDED M AND mu FOR l1-REGULARIZED WEIGHTED LEAST SQUARES
%
%--------------------------------------------------------------------------
% DESCRIPTION AND INPUTS
%--------------------------------------------------------------------------
% 
% [M,mu,A,b,sig,kap,tau,M12] = getM_mu(full,mu,m,n,Ameth,A,b,sig1,sig2,alpha)
%
% Constructs M and mu as described in Hale, Yin and Zhang 2007, based
% on the noise estimates sig1 and sig2, and the statistical parameter
% alpha.  
%
%  full  - only relevant if sig1 > 0 and AA' is not a multiple of I.  In 
%          this case, and if full == true, then the full estimate for M, 
%          M = (sig1^2 AA' + sig2^2 I)^{-1} is calculated.  Otherwise M is 
%          estimated to be (sig1^2 lam_max(AA') + sig2^2)^{-1} I.  The 
%          constant is then moved to mu so that M is returned as [].
%
%  mu    - if empty, then recommended mu value is calculated.
%
%  m,n   - A is an m x n matrix.
%
%  Ameth - indicates what type of A matrix is being passed.  values >= 0
%          correspond to indices from getData.m.
%               -1 - generic A matrix
%                0 - randn(m,n) (uses analytic bound on min and max
%                    eigenvalues of AA')
%                1 - 0 and then columns are scaled to unit norm
%                2 - 0 and then QR factorization to obtain orthonormal rows
%                    (AA' = I)
%                3 - bernoulli +/- 1 distribution
%                4 - partial Hadamard matrix (AA' = nI)
%                5 - partial fourier matrix (implicit,AA'=nI,ifft = A'/n)
%                6 - partial discrete cosine matrix (implicit,AA'=I)
%                7 - partial 2-d discrete cosine matrix (as 6)
%          Codes -1, 1, and 3 are equivalent; 2, 6 and 7 are equivalent; 4
%          and 5 are equivalent.
%
%  A,b   - problem data.  Only accessed if Ameth <= 4 (A is explicit).
%
%  sig1, sig2 - noise level estimates.  In particular we assume that 
%          b = A*(xs + epsilon1) + epsilon2, where epsiloni is iid Gaussian
%          noise with standard deviation sigi.
%
%  alpha - mu estimate uses 1-alpha chi^2 critial value.  mu is not very
%          sensitive to alpha and alternative is provided in case the 
%          statistics toolbox is not available.  Results in Hale, Yin and 
%          Zhang 2007 used alpha = 0.5.
%
%--------------------------------------------------------------------------
% OUTPUTS
%--------------------------------------------------------------------------
%
% M     - recommended value of M.  Returned as [] if M = I.
%
% mu    - if input mu is the empty matrix, output mu is the recommended 
%         value.  otherwise, mu is untouched.
%
% A,b   - problem data.  Untouched if AA' = I or A is implicit.  Otherwise, 
%         A and b are scaled so that the maximum eigenvalue of A'*M*A = 1.
%
% sig   - approximate noise level (standard deviation) of least squares 
%         solution.  equal to sqrt(sig1^2 + sig2^2/(min eig of AA')).
%
% kap   - condition number of M12*AA'*M12.
%
% tau   - if AtMA is not well conditioned, tau is returned as 2-eps.  
%         otherwise, tau is returned empty to signal that fpc default 
%         should be used.
%
% M12   - equal to M^(1/2) if ~isempty(M).
%--------------------------------------------------------------------------

function [M,mu,A,b,sig,kap,tau,M12] = getM_mu(full,mu,m,n,Ameth,A,b,sig1,sig2,alpha)

tau = []; M12 = [];

% options for eigs
opts.tol = 1E-4;
opts.disp = 0;
opts.issym = true;

% calculate M, mu, sig, kap, tau, M12
if full && sig1 > 0 && ismember(Ameth,[-1 0 1 3])
    % most general case--must calculate full M matrix
    AAt = A*A';
    M = sig1^2*AAt; M(1:m+1:m^2) = M(1:m+1:m^2) + sig2^2;
    % invert to get M
    lsopts.SYM = true; lsopts.POSDEF = true;
    M = linsolve(M,eye(m),lsopts);
    % need matrix square root
    M12 = sqrtm(M);
    M12AAtM12 = M12*AAt*M12;
    % get eigenvalues
    smax = sqrt(eigs(M12AAtM12,1,'lm',opts));
    smin = sqrt(eigs(M12AAtM12,1,0,opts));
    if Ameth == 0
        sig = sqrt(sig1^2 + sig2^2/max((1 - sqrt(m/n))^2*n,eps^2));
    else
        sig = sqrt(sig1^2 + sig2^2/max(eigs(AAt,1,0,opts),eps^2));
    end
    kap = smax^2/smin^2;
    % no tau because full M means well-conditioned AtMA
    if isempty(mu)
        if strcmp(which('chi2inv'),'');
            % avoids chi2inv, is good approximation for large n
            mu = smax*sqrt(n*kap/m);
        else
            mu = smax*sqrt(n*kap/chi2inv(1-alpha,m));
        end
    end
else
    % M is a multiple of I
    M = [];
    switch Ameth
        case {-1,1,3}
            % unknown eigenvalues
            AAt = A*A';
            smax = sqrt(eigs(AAt,1,'lm',opts));
            smin = sqrt(eigs(AAt,1,0,opts));
            muSig = sqrt(sig1^2*smax^2 + sig2^2);
            sig = sqrt(sig1^2 + sig2^2/smin^2);
            kap = smax^2/smin^2;
            tau = 2-eps;
        case 0
            % have tight bounds
            delta = m/n;
            smax = (1 + sqrt(delta))*sqrt(n);
            smin = max((1 - sqrt(delta))*sqrt(n),eps);
            muSig = sqrt(sig1^2*smax^2 + sig2^2);
            sig = sqrt(sig1^2 + sig2^2/smin^2);
            kap = smax^2/smin^2;
            tau = 2-eps;
        case {2,6,7}
            % AAt = I
            smin = 1; smax = 1;
            sig = sqrt(sig1^2 + sig2^2);
            muSig = sig;
            kap = 1;
        case {4,5}
            % AAt = nI
            smin = sqrt(n); smax = sqrt(n); 
            muSig = sqrt(sig1^2*n + sig2^2);
            sig = sqrt(sig1^2 + sig2^2/n);
            kap = 1;
            % if Ameth = 5, b was already scaled by getData since A scaling
            % is implemented in pfft_wrap_part
    end
    % calculate recommended mu if required
    if isempty(mu)
        if strcmp(which('chi2inv'),'');
            % avoids chi2inv, is good approximation for large n
            mu = (smax/muSig)*sqrt(n*kap/m);
        else
            mu = (smax/muSig)*sqrt(n*kap/chi2inv(1-alpha,m));
        end
    end
end

% scale A and b if required
if ismember(Ameth,[-1,0,1,3,4])
    A = (1/smax)*A; b = (1/smax)*b;
end

return

% Copyright (c) 2007.  Elaine Hale, Wotao Yin, and Yin Zhang
%
% Last modified 28 August 2007.