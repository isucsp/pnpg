classdef ADMM_NNL1 < Methods
    properties
        stepShrnk = 0.5;
        maxItr=1e2;
        s
        Psi_s
        Psi_sy
        pa
        rho;
        y;
        r_norm
        absTol=1e-4;
        preS
        cnt
        adaptiveStep=true;
        cumu=0;
        cumuTol=4;
    end
    methods
        function obj = ADMM_NNL1(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            alpha(alpha<0)=0;
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi;
            obj.Psit = Psit;
            obj.s = obj.Psit(alpha);
            obj.Psi_s = alpha;
            obj.coef(1) = 1;
            fprintf('use ADMM_NNL1 method\n');

            obj.rho=1;
            obj.y=0;
            obj.preS=obj.s;
            obj.cnt=0;
        end
        function main(obj)
            obj.warned = false;
            pp=0; obj.debug='';

            while(pp<obj.maxItr)
                obj.p = obj.p+1;
                pp=pp+1;

                [oldCost,obj.grad] = obj.func(obj.alpha);

                temp=obj.rho*obj.Psi(obj.s+obj.y);

                % start of line Search
                obj.ppp=0; goodStep=true; incStep=false; goodMM=true;
                while(true)
                    if(obj.adaptiveStep && ~incStep && obj.cumu>=obj.cumuTol)
                        % adaptively increase the step size
                        obj.t=obj.t*obj.stepShrnk;
                        obj.cumu=0;
                        incStep=true;
                    end
                    obj.ppp = obj.ppp+1;

                    newX = (obj.t*obj.alpha-obj.grad+temp)/(obj.t+obj.rho);
                    newX(newX<0)=0;

                    [newCost,newGrad]=obj.func(newX);
                    
                    if(innerProd(newX-obj.alpha,newGrad-obj.grad)<=sqrNorm(newX-obj.alpha)*obj.t/2)
                        if(obj.p<=obj.preSteps && obj.ppp<18 && goodStep && obj.t>0)
                            obj.t=obj.t*obj.stepShrnk; continue;
                        else
                            break;
                        end
                    else
                        if(obj.ppp<=20 && obj.t>0)
                            obj.t=obj.t/obj.stepShrnk; goodStep=false; 
                            if(incStep)
                                obj.cumuTol=obj.cumuTol+4;
                                incStep=false;
                            end
                        else
                            goodMM=false;
                            obj.debug=[obj.debug 'falseMM'];
                            break;
                        end
                    end
                end
                obj.stepSize = 1/obj.t;
                obj.fVal(3) = obj.fArray{3}(newX);
                obj.cost = newCost+obj.u*obj.fVal(3);

                obj.difAlpha = relativeDif(obj.alpha,newX);
                obj.alpha = newX;

                if(obj.ppp==1 && obj.adaptiveStep)
                    obj.cumu=obj.cumu+1;
                else
                    obj.cumu=0;
                end

                PsitAlpha=obj.Psit(obj.alpha);

                obj.s = Utils.softThresh(PsitAlpha-obj.y,obj.u/obj.rho);
                obj.y = obj.y + (obj.s-PsitAlpha);
                
                difS=pNorm(obj.s-obj.preS); obj.preS=obj.s;
                residual = pNorm(obj.s-PsitAlpha);
                sNorm = pNorm(obj.s);

                if(difS<=obj.absTol*sNorm && residual<=obj.absTol*sNorm) break; end
                if(obj.cnt>10) % prevent back and forth adjusting
                    if(difS>10*residual)
                        obj.rho = obj.rho/2; obj.y=obj.y*2; obj.cnt=0;
                    elseif(difS<residual/10)
                        obj.rho = obj.rho*2; obj.y=obj.y/2; obj.cnt=0;
                    end
                end
            end
        end

    end
end

