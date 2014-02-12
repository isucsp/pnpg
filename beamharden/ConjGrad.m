classdef ConjGrad < handle
    properties
        cgType = 'prp';
        func
        fArray
        coef
        preG=1;
        preP=0;
        converged = false;
    end
    methods
        function obj = ConjGrad()
            fval = zeros(3,length(fArray));
            if(~isfield(opt,'t3'))
                [temp,temp1]=polyIout(0,Ie);
                t3=max(abs(PsitPhitz+PsitPhit1*log(temp)))*temp1/temp;
                t3=t3*10^opt.a; out.t3 = t3;
            end
        end
        function main(obj)
            pp=0; obj.converged = false; obj.warned = false;
            while(pp<obj.maxStepNum)
                pp=pp+1;
                [oldCost,grad,hessian] = obj.func(obj.alpha);

                if(min(weight)<=0)
                    warning(['fAlpha: obj function is non-convex over alpha,'...
                        '(%g,%g)'],min(weight),max(weight));
                    weight(weight<0)=0;
                    str='';
                end

                maxStep=1;
                beta=grad'*(grad-preG)/(preG'*preG);
                deltaAlpha=grad+max(beta,0)*preP;
                deltaNormAlpha=grad'*deltaAlpha;
                s1=deltaNormAlpha/h(deltaAlpha,2);
                %atHessianA(deltaAlpha,weight,t1*hphi,Phi,Phit,t3, Psit,opt.muLustig,sqrtSSqrMu);
                preP = deltaAlpha; preG = grad;
                deltaAlpha = deltaAlpha*s1;
                deltaNormAlpha = deltaNormAlpha*s1;

                if(deltaNormIe<obj.thresh)
                    obj.converged=true;
                end

                % start of line Search
                ppp=0; stepSz=min(1,maxStep);
                while(~obj.converged)
                    ppp=ppp+1;
                    newX=alpha-stepSz*deltaAlpha;
                    %newX(newX<0)=0; % force it be positive;
                    newCost=obj.func(newX);

                    if(newCost <= oldCost - stepSz/2*deltaNormAlpha)
                        obj.difAlpha(p) = norm(obj.alpha(:)-newX(:))^2;
                        obj.alpha = newX;
                        break;
                    else
                        if(ppp>10)
                            out.llAlphaDif(p) = 0;
                            warning('exit iterations for higher convergence criteria: %g\n',deltaNormIe);
                            obj.cost = oldCost;
                            obj.warned = true;
                            obj.converged = true;
                            
                        else
                            stepSz=stepSz*stepShrnk;
                        end
                    end
                end
                %end of line search
            end
        end
        function x= cg(c,hessianA,atHessianA,maxItr)
            % This function solve the problem 
            % min c'*x+1/2 atHessianA(x)
            % hessianA=hessian*x
            x=0; g=c; p=0; i=0;
            while(i<maxItr)
                i= i+1;
                preP=p; preG=g;
                g=c+hessianA(x);
                p=-g-g'*g/(preG'*preG)*preP;
                x=x-p*(p'*g/atHessianA(p));
            end
        end

        function oldJunk()
            if(t1==0 && p==1)
                if(interiorPointAlpha)
                    t1=min([1, abs(diff0'*difphi/norm(difphi)), abs(costA/costB)]);
                else t1=1;
                end
            end

            if(0)
                %afun=@(xxx,yyy) fhessianA(xxx);
                %[deltaAlpha,~]=bicg(afun,difAlpha,[],5);
                %afun=@(xxx,yyy) fhessianA(xxx);
                fhessianA=@(gAlpha) hessianA(gAlpha,weight,hphi*t1,Phi,Phit);
                fatHessianA=@(gAlpha) atHessianA(gAlpha,weight,t1*hphi,Phi,Phit);
                deltaAlpha = cg(-difAlpha,fhessianA,fatHessianA,1);
            end
            if(interiorPointAlpha)
                temp=find(deltaAlpha>0);
                if(isempty(temp)) maxStep=1;
                else maxStep=min(alpha(temp)./deltaAlpha(temp)); end
                maxStep=maxStep*0.99;
            end
            if(interiorPointAlpha)
                if(deltaNormAlpha<1e-5)
                    if(t1 < 1e-10/length(alpha))
                        alphaReady=1;
                    else
                        t1=t1/10;
                        thresh1=thresh1/10;
                    end
                end
            end

        end
    end

