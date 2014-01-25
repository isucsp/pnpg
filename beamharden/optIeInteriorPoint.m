f1=@(IeVar) barrierIe(IeVar);
[costC,diff1,h1]=f1(Ie);

if(t2==0)
    t2=min([1, abs(diff0'*diff1/norm(diff1)), costA/costC]);
    temp=ceil(log(t2/t2Lim*E));
    thresh2q=10; %(thresh2/thresh2Lim)^(1/temp);
end
%if(t2<=1e-5) maxPP=100; end
cost1=costA+t2*costC;
dif=diff0+t2*diff1;
H=h0+t2*h1;
temp=eig(H);
if(all(imag(temp)==0))
    if(min(temp)<=0)
        %save('negativeDefiniteHessian.mat');
        fprintf([repmat(' ',1,80) '\n']);
        fprintf(['p=%d, pp=%d, Heissein is not positive '...
            'definite, minnimum=%g\n'], p, pp, min(temp));
        fprintf([repmat(' ',1,80) '\n']); str='';
        temp=dif'*H*dif;
        if(temp>0) deltaIe=dif'*dif/temp*dif;
        else
            fprintf('after correction, still negative, give up\n');
            out='failed';
            return;
        end
    else
        deltaIe=H\dif;
    end
else
    save('noRealEigen.mat');
    fprintf([repmat(' ',1,80) repmat('\n',1,1)]);
    fprintf('Error for the H matrix!!\n');
    fprintf('p=%d, pp=%d, Heissein has complex eigen value\n',...
        p, pp);
    fprintf(repmat([repmat(' ',1,80) repmat('\n',1,1)],1,3)); 
    str='';
    return;
end

temp=find(deltaIe>0);
if(isempty(temp)) opt.maxStep=1.1;
else opt.maxStep=min(Ie(temp)./deltaIe(temp)); end

if(sum(deltaIe)<0)
    opt.maxStep=min(opt.maxStep, (sum(Ie)-1)/sum(deltaIe));
end
opt.maxStep=opt.maxStep*0.99;

opt.g=dif;
penalty=@(x) t2*f1(x);
out = lineSearch(Ie,deltaIe,f0,penalty,1,cost1,opt);
res2(p)=out.costA; penaltyIe=out.costB; stepSzIe=out.stepSz;
deltaNormIe=out.delta; Ie=out.x; zmf=out.zmf;
%if(out.stepSz<1e-10) break; end
if(deltaNormIe<thresh2)
    if(t2 < t2Lim/(E+1))
        IeReady=1;
        break;
    else
        t2=t2/10;
        thresh2=thresh2/thresh2q;
    end
end
