function Psi=wvltMat(p,dwt_L,daub)

    if(nargin<3)
        daub = 4;
        if(nargin<2)
            dwt_L= 3;        %levels of wavelet transform
        end
    end
    wvltName = sprintf('PsiWvlt%dL%dD%d.mat',p,dwt_L,daub);
    if(exist(wvltName,'file'))
        load(wvltName);
        return;
    end
    % caution: don't use the wavelet tools from matlab, it is slow
    daub_H = daubcqf(daub);
    psi = @(xx) midwt(xx,daub_H,dwt_L);
    Psi = zeros(p,p);
    x=zeros(p,1);
    str='';
    for j=1:p
        strlen=length(str); str=sprintf('j=%d',j); fprintf([repmat('\b',1,strlen) '%s'],str);
        x(j)=1;
        Psi(:,j)=psi(x);
        x(j)=0;
    end

    save(wvltName,'Psi');
end
