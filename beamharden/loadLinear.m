function opt=loadLinear(obj,opt)
    %s = RandStream.create('mt19937ar','seed',0);
    %RandStream.setGlobalStream(s);

    tics=[0,...
        90, 190,...
        210, 240,...
        260, 320, 340, 370,...
        440, 490, 540,...
        610, 640, 650,...
        700, 750, 780,...
        820, 880, 940,...
        1023];
    t=(0:1023); t=t/max(t)*max(tics);
    x(t>=tics( 1) & t<tics( 2))=0;
    x(t>=tics( 2) & t<tics( 3))=linspace(0,1.4,sum(t>=tics(2) & t<tics(3)));
    x(t>=tics( 3) & t<tics( 4))=0;
    x(t>=tics( 4) & t<tics( 5))=1.5*0;
    x(t>=tics( 5) & t<tics( 6))=0;
    x(t>=tics( 6) & t<tics( 7))=2.3;
    x(t>=((tics( 6)*2+tics(7))/3) & t<tics( 7))=linspace(2.3,1.6,sum(t>=((tics( 6)*2+tics(7))/3) & t<tics( 7)));
    x(t>=tics( 7) & t<tics( 8))=0.3;
    x(t>=tics( 8) & t<tics( 9))=1.0;
    x(t>=tics( 9) & t<tics(10))=0;
    x(t>=tics(10) & t<tics(11))=linspace(0.9,1.7,sum(t>=tics(10) & t<tics(11)));
    x(t>=tics(11) & t<tics(12))=linspace(1.7,0.9,sum(t>=tics(11) & t<tics(12)));
    x(t>=tics(12) & t<tics(13))=0;
    x(t>=tics(13) & t<tics(14))=linspace(1.1,1.1,sum(t>=tics(13) & t<tics(14)));
    x(t>=tics(14) & t<tics(15))=1.5;
    x(t>=tics(15) & t<tics(16))=0; %1.5+0.6*sin((0:sum(t>=tics(15) & t<tics(16))-1)/sum(t>=tics(15) & t<tics(16))*pi);
    x(t>=tics(16) & t<tics(17))=1.5;
    x(t>=tics(17) & t<tics(18))=linspace(1.5,1.1,sum(t>=tics(17) & t<tics(18)));
    x(t>=tics(16) & t<tics(18))=sin((1:sum(t>=tics(16) & t<tics(18)))/sum(t>=tics(16) & t<tics(18))*pi);
    x(t>=tics(18) & t<tics(19))=0;
    x(t>=tics(19) & t<tics(20))=1.6*(((1:sum(t>=tics(19) & t<tics(20)))/sum(t>=tics(19) & t<tics(20))).^2);
    x(t>=tics(20) & t<tics(21))=1.6*(((sum(t>=tics(20) & t<tics(21))-1:-1:0)/sum(t>=tics(20) & t<tics(21))).^2);
    x(t>=tics(21))=0;

    if(~exist('obj'))
        [~,l1norm]=Utils.testSparsity(x);
        opt.l1norm=l1norm;
        figure; plot(x);
        return;
    end

    %figure(1); plot(t,x); ylim([-2,5]);
    if(~exist('opt')) opt.m=500; end
    if(~isfield(opt,'m')) opt.m=500; end
    if(~isfield(opt,'snr')) opt.snr=inf; end
    if(~isfield(opt,'noiseType')) opt.noiseType='gaussian'; end
    if(~isfield(opt,'matrixType')) opt.matrixType='gaussian'; end
    if(~isfield(opt,'padZero')) opt.padZero=0; end
    x=[x(:); zeros(opt.padZero,1)];
    n = length(x);      % number of features

    x0=x(:);
    if(strcmpi(opt.matrixType,'nonneg'))
        A = rand(opt.m,n);
        a=0.3;
        idx=(A<a);
        A(idx)=0;
        A(~idx)=(A(~idx)-a)/(1-a);
    elseif(strcmpi(opt.matrixType,'gaussian'))
        A = randn(opt.m,n);
        A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
    else
        error('error input matrixType');
    end

    if(strcmpi(opt.noiseType,'gaussian'))
        v = randn(opt.m,1);
        v = v*(norm(A*x0)/sqrt(opt.snr*opt.m));
        b = A*x0 + v;
        % fprintf('nnz(x0) = %d; signal-to-noise ratio: %.2f\n', nnz(x0), norm(A*x0)^2/norm(v)^2);
    elseif(strcmpi(opt.noiseType,'poisson'))
        y = A*x0;
        a = (opt.snr-1)*sum(y)/sum(y.^2);
        a = 1000; % this is equivalent to SNR ~ 3e5
        %x0=x0*1000;
        x0=a*x0;
        b = poissrnd(a*y);
        % figure; showImg(A);
    else
        error('wrong input noiseType');
    end

    %fprintf('solving instance with %d examples, %d variables\n', opt.m, n);

    gamma_max = norm(A'*b,'inf');
    gamma = 0.1*gamma_max;

    obj.trueImg = x0;
    obj.y = b;

    obj.Phi = @(xx) A*xx(:);
    obj.Phit = @(xx) A'*xx(:);

    obj.dwt_L= 3;        %levels of wavelet transform
    obj.daub = 4;

    % caution: don't use the wavelet tools from matlab, it is slow
    daub_H = daubcqf(obj.daub);
    obj.Psi = @(xx) midwt(xx,daub_H,obj.dwt_L);
    obj.Psit= @(xx)  mdwt(xx,daub_H,obj.dwt_L);

    opt.trueAlpha=x0;
    opt.L = max(svd(A))^2;
    if(strcmpi(opt.noiseType,'poisson'))
        opt.L=opt.L/min(obj.y);
    end
end

