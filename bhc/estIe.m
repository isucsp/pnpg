function Ie = estIe(kappa,Ie,)
    currNumIe=length(kappa);
    for numIe=currNumIe:-1:1
        temp=logspace(log10(kappa(1)),log10(kappa(end)),numIe);
        Ie=exp(interp1(log(massAttenCoef(idx1:idx2,1)),...
            log(massAttenCoef(idx1:idx2,2)), log(epsilon),'spline'));

        kappa=temp(:);  %*mean(X(find(idx(:)==i+1))); %/(1-(opt.K-1)*eps);
        q=kappa(2)/kappa(1);
        temp = [temp(1)^2/temp(2); temp(:); temp(end)^2/temp(end-1)];
        polymodel = Spline(opt.spectBasis,temp);
        polymodel.setPlot(opt.trueKappa,opt.trueIota,opt.epsilon);
        polyIout = polymodel.polyIout;

        Ie=Ie/polyIout(0,Ie);

        IeStep = NPG_ind(Ie,true,[],[],opt.maxIeSteps,opt.stepShrnk,opt.thresh);
        A = polyIout(Phi(alpha),[]);
        IeStep.func = @(III) gaussLI(Imea,A,III);
        if(IeStep.Ie(1) || IeStep.Ie(end))
            fprintf('\nExtend the spectrum to the');
            if(IeStep.Ie(1))
                kappa=[kappa(1)/q; kappa];
                IeStep.setIe([0; IeStep.Ie]);
                Ie=[0; Ie];
                fprintf(' left');
            end
            if(IeStep.Ie(end))
                kappa=[kappa; kappa(end)*q];
                IeStep.setIe([IeStep.Ie; 0]);
                Ie=[Ie; 0];
                fprintf(' right');
            end
            temp = [kappa(1)^2/kappa(2); kappa(:); kappa(end)^2/kappa(end-1)];
            polymodel = Spline(opt.spectBasis,temp);
            polymodel.setPlot(opt.trueKappa,opt.trueIota,opt.epsilon);
            polyIout = polymodel.polyIout;
            fprintf(' !\n'); strlen=0;
        end

        IeStep.main();
    end
