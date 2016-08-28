function u=uBound(Psi,Psit,Pncx,xstar,g)
    % 
    % g is the gradient of L(x) at xstar
    % author: Renliang Gu
    %

    g=g(:);
    xstar=xstar(:);
    v=0;
    t=zeros(size(xstar));
    w=zeros(size(Psit(xstar)));
    z=zeros(size(t));
    Psi_w=Psi(w);

    rho=1;

    cnt=0;
    ii=0;

    opts.printEvery=inf;
    opts.maxTotalIts=5e3;
    opts.maxIts=1000;
    opts.factr=1e-6/eps;

    normG=norm(g,'fro');
    relDif=1;

    lb=-ones(size(w)); ub=ones(size(w));
    ww=sign(Psit(xstar)); A=(ww~=0);
    lb(A)=ww(A); ub(A)=ww(A);

    while(relDif>=1e-9)

        ii=ii+1; cnt=cnt+1;

        preV=v; preT=t; preW=w; preZ=z; prePsi_w=Psi_w;
         
        % objective: 0.5*||Psi(w)+t+z+v*g||_F^2
        % subject to the [-1,1] box
        opts.x0=w;
        [w,cost,info]=lbfgsb(@(x) quadbox(x,v*g+t+z,Psi,Psit),...
            lb,ub,opts);
        numItr=info.iterations;

        Psi_w=Psi(w);


        v=(1-rho*g'*(Psi_w+t+z)) / (rho*(g'*g));
        t=Pncx( -v*g-Psi_w-z );

        z=z+Psi_w+v*g+t;

        pRes=norm(v*g+Psi_w+t,2);
        dRes1=g'*(prePsi_w-Psi_w+preT-t);
        dRes2=norm(prePsi_w-Psi_w,2);
        gap=z'*(v*g+Psi_w+t);

        relDif=max(abs([pRes, dRes1, dRes2]))/normG;

        fprintf('itr=%d, u=%g relDif=%g gap=%g\n',ii, 1/v, relDif, gap);
        %fprintf('itr=%d, u=%g %g residual=%g rho=%g difPQ=%g, gap=%g numItr=%d normG=%g\n',...
        %   ii, max(abs([p(:); q(:)])), u, residual, rho, difPQ,gap, numItr, normG);

        %if(cnt>10) % prevent excessive back and forth adjusting
        %    if(difPQ>10*residual)
        %        rho=rho/2 ; cnt=0;
        %    elseif(difPQ<residual/10)
        %        rho=rho*2 ; cnt=0;
        %    end
        %end
    end
    u=1/v;

    function [f,g] = quadbox(x,y,Psi,Psit)
        r=Psi(x)+y;
        f=0.5*norm(r,2)^2;
        if(nargout>1)
            g=Psit(r);
        end
    end
end
