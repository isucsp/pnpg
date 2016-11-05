function u=uBound(Psi,Psit,Pncx,xstar,g)
    % 
    % g is the gradient of L(x) at xstar
    % author: Renliang Gu
    %

    g=g(:);

    if(norm(g+Pncx(-g))==0) 
        u=0;
        return;
    end

    % it is sure that normG!=0
    normG=norm(g,2);
    g=g/normG;

    xstar=xstar(:);
    v=1;
    vg=v*g;
    t=zeros(size(xstar));
    w=zeros(size(Psit(xstar)));
    z=zeros(size(t));
    Psi_w=Psi(w);

    rho=1;

    cnt=0;
    ii=0;
    EndCnt=0;

    dRes=1e-4;

    opts.printEvery=inf;
    opts.maxTotalIts=5e3;
    opts.maxIts=1000;
    opts.factr=0; %1e-6/eps;

    strlen=0;

    lb=-ones(size(w)); ub=ones(size(w));
    ww=sign(Psit(xstar)); A=(ww~=0);
    lb(A)=ww(A); ub(A)=ww(A);

    while(EndCnt<3 && ii<=5e3)

        ii=ii+1; cnt=cnt+1;

        preV=v; preT=t; preW=w; preZ=z; prePsi_w=Psi_w;
         
        opts.pgtol=max(1e-2*dRes,1e-14);
        for j=1:1

            % objective: 0.5*||Psi(w)+t+z+v*g||_F^2
            % subject to the [-1,1] box
            opts.x0=w;
            [w,cost,info]=lbfgsb(@(x) quadbox(x,vg+t+z,Psi,Psit),...
                lb,ub,opts);
            %numItr=info.iterations; %disp(info); strlen=0;
            Psi_w=Psi(w);

            % Since ||g||=1, v=(1-rho*g'*(Psi_w+t+z)) / (rho*(g'*g))
            % reduces to
            v=1/rho-g'*(Psi_w+t+z); vg=v*g;
            t=Pncx( -vg-Psi_w-z );
        end

        z=z+Psi_w+vg+t;

        pRes=norm(vg+Psi_w+t,2);
        %dRes1=g'*(prePsi_w-Psi_w+preT-t);
        %dRes2=norm(prePsi_w-Psi_w,2)^2;
        dRes1=abs(preV-v);  %norm(g*(preV-v));
        dRes2=norm(preT-t);
        dRes=[max(dRes1,dRes2)];
        gap=z'*(vg+Psi_w+t);

        str=sprintf('itr=%d, u=%g pRes=%g dRes1=%g dRes2=%g gap=%g rho=%g                 ',ii, normG/v, pRes,...
                dRes1, dRes2, gap, rho);
        if(strlen==0 || (mod(ii-1,100)==0 || (ii<=100 && mod(ii-1,10)==0) || ii-1<10))
            fprintf('\n%s',str);
        else
            fprintf([repmat('\b',1,strlen) '%s'],str);
        end
        strlen = length(str);

        if(cnt>100) % prevent excessive back and forth adjusting
            if(dRes>10*pRes)
                rho=rho/2; z=z*2; cnt=0;
            elseif(dRes<pRes/10)
                rho=rho*2; z=z/2; cnt=0;
            end
        end

        if(max(pRes,dRes)<1e-9)
            EndCnt=EndCnt+1;
        else
            EndCnt=0;
        end
    end
    u=normG/v;
    fprintf('\n\n');

    function [f,g] = quadbox(x,y,Psi,Psit)
        r=Psi(x)+y;
        f=0.5*norm(r,2)^2;
        if(nargout>1)
            g=Psit(r);
        end
    end
end

