function nestroAnim(func)

    if(~exist('func','var')) func='Poisson'; end
    switch lower(func)
        case lower('poisson')
            A=[10,0.21;2000.5,5000];
            x0=[0.4;0.4];
            y=poisson(A*x0*100)/100;
            y=A*x0;
            f=@(a,b) (A(1,1)*a+A(1,2)*b-y(1)+A(2,1)*a+A(2,2)*b-y(2))...
                +y(1)*log(y(1)./(A(1,1)*a+A(1,2)*b))...
                +y(2)*log(y(2)./(A(2,1)*a+A(2,2)*b));
            [x1,x2]=meshgrid(linspace(0,3*x0(1),1000), linspace(0,3*x0(2),1000));
            fx=f(x1,x2);
            Phi=@(x) A*x;
            Phit=@(x) A.'*x;

            %xInit=Phit(y);
            xInit=[0.1,1];
            xInit=[6,5];
            xInit=[0,3];
            xInit=[1,0];
            xInit=[5.5,0.5];
            xInit=[4,5];
            xInit=[3,0];
            xInit=[1,1]*1e-1;

            opt.noiseType='Poisson';

            v=logspace(-2,4,20);

            opt.L = max(svd(A))^2;
        case lower('rosenbrock')
            opt.noiseType='other';
            opt.f=@(varargin) rosenbrock(1,100,varargin{:});

            x0=[1.2;1];
            [x1,x2]=meshgrid(linspace(0,3*x0(1),1000), linspace(0,3*x0(2),1000));
            fx=opt.f(x1,x2);
            v=logspace(log10(max(1e-3,min(fx(:)))),log10(max(fx(:))),20);
            Phi=[]; Phit=[]; y=[];
            xInit=[1; 3];
    end

    Psi=@(x) x;
    opt.maxItr=1e4; opt.thresh=1e-5; opt.debugLevel=1; opt.u=1e-10;
    opt.continuation=false; opt.saveXtrace=true;

    opt.alphaStep='pg'; opt.adaptiveStep=false;  opt.L=opt.L/10;
    out=solver(Phi,Phit,Psi,Psi,y,xInit,opt);
    out.alpha
    figure; contour(x1,x2,fx,v); hold on; plot(xInit(1),xInit(2),'r.');
    plot(out.alphaTrace(1,:),out.alphaTrace(2,:),'g');

    keyboard

    opt.alphaStep='pnpg'; opt.adaptiveStep=true;
    out=solver(Phi,Phit,Psi,Psi,y,xInit,opt);
    out.alpha
    figure; contour(x1,x2,fx,v); hold on; plot(xInit(1),xInit(2),'r.');
    plot(out.alphaTrace(1,:),out.alphaTrace(2,:),'r');

    keyboard

    opt.alphaStep='pnpg'; opt.adaptiveStep=false;
    out=solver(Phi,Phit,Psi,Psi,y,xInit,opt);
    out.alpha
    figure; contour(x1,x2,fx,v); hold on; plot(xInit(1),xInit(2),'r.');
    plot(out.alphaTrace(1,:),out.alphaTrace(2,:),'b');

end

function [f,g,h] = rosenbrock(a,b,x,y)
    % The Rosenbrock function
    % f(x)=(a-x(1))^2+b*(x(2)-x(1)^2)^2
    % g(x)_1=4b*x(1)^3-4*b*x(2)*x(1)+2*(x(1)-a)
    % g(x)_2=2*b*(x(2)-x(1)^2)

    c=2;
    if(nargin==4)
        f=(a-(x-c)).^2+b*(y-(x-c).^2).^2;
    else
        f=(a-x(1,:)+c).^2+b*(x(2,:)-(x(1,:)-c).^2).^2;
    end

    if(nargout>=2)
        g=zeros(size(x));
        g(1,:)=4*b*(x(1,:)-c).^3-4*b*x(2,:).*(x(1,:)-c)+2*(x(1,:)-c-a);
        g(2,:)=2*b*(x(2,:)-(x(1,:)-c).^2);

        if(nargout>=3)
            H=[12*b*(x(1)-c)^2-4*b*x(2)+2, -4*b*(x(1)-c);
            -4*b*(x(1)-c), 2*b];
            h=@(x,opt) hessian(H,x,opt);
        end
    end

    function hh=hessian(H,x,opt)
        if(opt==1)
            hh = H*x;
        else
            hh = x'*H*x;
        end
    end
end




