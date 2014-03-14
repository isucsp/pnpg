function opt=loadLasso(obj,opt)
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);

    m = 500;       % number of examples
    n = 2500;      % number of features

    x0 = sprandn(n,1,0.05);
    x0 = abs(x0);
    A = randn(m,n);
    A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
    v = sqrt(0.001)*randn(m,1);
    b = A*x0 + v;

    fprintf('solving instance with %d examples, %d variables\n', m, n);
    fprintf('nnz(x0) = %d; signal-to-noise ratio: %.2f\n', nnz(x0), norm(A*x0)^2/norm(v)^2);

    gamma_max = norm(A'*b,'inf');
    gamma = 0.1*gamma_max;

    % cached computations for all methods
    AtA = A'*A;
    Atb = A'*b;

    obj.trueImg = x0;
    obj.mask = 1:length(x0);
    obj.mask = obj.mask(:);
    obj.y = b;

    obj.Phi = @(xx) A*xx(:);
    obj.Phit = @(xx) A'*xx(:);
    obj.FBP = @(xx) A'*xx(:);
    obj.Psi = @(xx) xx(:);
    obj.Psit = @(xx) xx(:);

    opt.u = gamma;
    opt.trueAlpha=x0;
end
