function opt=loadLinear(obj,opt)
    %s = RandStream.create('mt19937ar','seed',0);
    %RandStream.setGlobalStream(s);

    tics=[0, 65, 190, 210, 240, 260, 320, 340, 370, 435, 490, 545, 610, 640, 650, 740, 750, 780, 820, 880, 940, 1023];
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
    x(t>=tics(13) & t<tics(14))=linspace(1.1,1.5,sum(t>=tics(13) & t<tics(14)));
    x(t>=tics(14) & t<tics(15))=1.5;
    x(t>=tics(15) & t<tics(16))=1.5+0.6*sin((0:sum(t>=tics(15) & t<tics(16))-1)/sum(t>=tics(15) & t<tics(16))*pi);
    x(t>=tics(16) & t<tics(17))=1.5;
    x(t>=tics(17) & t<tics(18))=linspace(1.5,1.1,sum(t>=tics(17) & t<tics(18)));
    x(t>=tics(18) & t<tics(19))=0;
    x(t>=tics(19) & t<tics(20))=1.8*((1:sum(t>=tics(19) & t<tics(20)))/sum(t>=tics(19) & t<tics(20))).^2;
    x(t>=tics(20) & t<tics(21))=1.8*((sum(t>=tics(20) & t<tics(21)):-1:1)/sum(t>=tics(20) & t<tics(21))).^2;
    x(t>=tics(21))=0;

    if(~exist('obj'))
        for i=0.7:0.01:0.98
            [~,l1norm]=test(x,i);
        end
        [i,j]=find(l1norm==min(l1norm(l1norm>0)));
        fprintf('level=%d,daub=%d  :  |x|_1 of signal achieves %g\n',i,2*j,l1norm(i,j));
        opt.l1norm=l1norm;
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
        x0=x0*1000;
        b = poissrnd(A*x0);
        % b=A*x0;
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

    obj.dwt_L= 5;        %levels of wavelet transform
    obj.daub = 4;

    % caution: don't use the wavelet tools from matlab, it is slow
    daub_H = daubcqf(obj.daub);
    obj.Psi = @(xx) midwt(xx,daub_H,obj.dwt_L);
    obj.Psit= @(xx)  mdwt(xx,daub_H,obj.dwt_L);

    opt.trueAlpha=x0;
end

function [dif,l1norm]=test(x,perce)
    for i=1:9
        for j=1:4
            level = i; wave = 2*j;
            [dif(i,j),~,l1norm(i,j)]=getError(x,perce,level,wave);
        end
    end
    [i,j]=find(dif==min(dif(dif>0)));
    fprintf('level=%d,daub=%d  :  %g%% of signal achieves %g%%\n',i,2*j,(1-perce)*100,dif(i,j)*100);
end

function [dif,coef,l1norm]=getError(x,perce,level,wave)
    h=daubcqf(wave);
    [wav,L] = mdwt(x,h,level);
    [~,idx]=sort(abs(wav));
    coef=wav(idx);
    coef=coef(end:-1:1);
    %figure(2); plot(wav(idx)); hold on;
    wav(idx(1:floor(perce*length(x))))=0;
    recX=midwt(wav,h,level);
    %figure(1); plot(recX); hold on; plot(x,'r');
    dif=(norm(x-recX)/norm(x))^2;
    l1norm=sum(abs(coef));
end

