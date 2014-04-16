function opt=loadLinear(obj,opt)
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);

    tics=[0, 40, 230, 254, 290, 310, 370, 390, 420, 470, 540, 610, 660, 690, 700, 790, 800, 830, 850, 920, 990, 1023];
    t=(0:1023); t=t/max(t)*max(tics);
    x(t>=tics( 1) & t<tics( 2))=0;
    x(t>=tics( 2) & t<tics( 3))=linspace(0,3,sum(t>=tics(2) & t<tics(3)));
    x(t>=tics( 3) & t<tics( 4))=0;
    x(t>=tics( 4) & t<tics( 5))=1.5;
    x(t>=tics( 5) & t<tics( 6))=0;
    x(t>=tics( 6) & t<tics( 7))=2.3;
    x(t>=((tics( 6)*2+tics(7))/3) & t<tics( 7))=linspace(2.3,1.6,sum(t>=((tics( 6)*2+tics(7))/3) & t<tics( 7)));
    x(t>=tics( 7) & t<tics( 8))=0.3;
    x(t>=tics( 8) & t<tics( 9))=1.0;
    x(t>=tics( 9) & t<tics(10))=0;
    x(t>=tics(10) & t<tics(11))=linspace(0.9,1.7,sum(t>=tics(10) & t<tics(11)));
    x(t>=tics(11) & t<tics(12))=linspace(1.7,0.9,sum(t>=tics(11) & t<tics(12)));
    x(t>=tics(12) & t<tics(13))=0;
    x(t>=tics(13) & t<tics(14))=linspace(1.1,1.5,sum(t>=tics(13) & t<tics(14)));
    x(t>=tics(14) & t<tics(15))=1.5;
    x(t>=tics(15) & t<tics(16))=1.5+0.6*sin((0:sum(t>=tics(15) & t<tics(16))-1)/sum(t>=tics(15) & t<tics(16))*pi);
    x(t>=tics(16) & t<tics(17))=1.5;
    x(t>=tics(17) & t<tics(18))=linspace(1.5,1.1,sum(t>=tics(17) & t<tics(18)));
    x(t>=tics(18) & t<tics(19))=0;
    x(t>=tics(19) & t<tics(20))=2.5*((1:sum(t>=tics(19) & t<tics(20)))/sum(t>=tics(19) & t<tics(20))).^2;
    x(t>=tics(20) & t<tics(21))=2.5*((sum(t>=tics(20) & t<tics(21)):-1:1)/sum(t>=tics(20) & t<tics(21))).^2;
    x(t>=tics(21))=0;

    %figure(1); plot(t,x); ylim([-2,5]);
    m = 500;       % number of examples
    n = length(x);      % number of features

    x0=x(:);
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

