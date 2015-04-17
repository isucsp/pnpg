function [CTdata, args] = genBeamHarden(symbol, densityMap, ops, varargin)
%   Simulate and generate measurements for beamhardening effect
%   Parameters includs:
%   'trueImg'
%   'iota'
%   'showImg'
%   'voltage'
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Fri 17 Apr 2015 03:52:46 PM CDT

    args = parseInputs(varargin{:});
    epsilon = args.epsilon(:);
    K=length(symbol);
    
    loadXrayMassCoef
    for k=1:K
        for i=1:length(symbols)
            if(strcmpi(symbol{k},symbols{i}))
                massAttenCoef=mac{i};
                massAttenCoef(:,1)=massAttenCoef(:,1)*1e3;

                density(k) = zaid(i,3);

                % from the original to 130 points, increase 19%
                % from 130 points to 1300 points, increasing 1e-6, which is not necessary.
                idx = massAttenCoef(:,1)>=min(epsilon) & massAttenCoef(:,1)<=max(epsilon);
                idx1= find(idx, 1 )-1; idx2 = find(idx, 1, 'last' )+1;
                kappa{k}=exp(interp1(log(massAttenCoef(idx1:idx2,1)),...
                    log(massAttenCoef(idx1:idx2,2)), log(epsilon),'spline'));

                if(args.showImg)
                    %figure; loglog(massAttenCoef(:,1),massAttenCoef(:,2));
                    %figure; semilogy(epsilon,kappa{k},'-');
                    %hold on;loglog(massAttenCoef(:,1),massAttenCoef(:,2),'*');
                    figure; loglog(massAttenCoef(idx1:idx2,2),...
                        gampdf((massAttenCoef(idx1:idx2,1)-20)*16/100,5,1),'*');
                    hold on; semilogy(kappa{k},args.iota);
                    %figure; semilogx(massAttenCoef(:,1),args.iota);
                    figure; semilogy(epsilon,args.iota,'*-');
                end

                break;
            end
        end
    end

    temp = [epsilon(1);(epsilon(1:end-1)+epsilon(2:end))/2;epsilon(end)];
    deltaEpsilon = temp(2:end)-temp(1:end-1);

    for k=1:K
        y{k}=ops.Phi(densityMap{k})*density(k);
    end

    Ts = 1;
    while(1)
        Imea=0;
        for i=1:length(epsilon)
            PhiAlpha = 0; 
            for k=1:K
                PhiAlpha = PhiAlpha + y{k} * kappa{k}(i); 
            end
            PhiAlpha=PhiAlpha*Ts;
            Imea=Imea+exp(-PhiAlpha)*args.iota(i)*deltaEpsilon(i);
        end

        temp = max(Imea(:));
        Imea = Imea/temp;
        if(min(Imea(:))<args.scale) Ts = Ts*0.8;
        else
            args.iota = args.iota/temp;
            break;
        end
    end

    fprintf('Sampling rate is determined to be %g\n',Ts);

    args.Ts = Ts;
    args.kappa = kappa{density==max(density)};
    args.density=max(density);
    CTdata=Imea;

    if(args.showImg)
        y=-log(Imea);
        rec=ops.FBP(y);
        rec2=ops.FBP(ops.Phi(densityMap{1}));
        figure; showImg(rec); figure; showImg(rec2);
    end
end

function args = parseInputs(varargin)

    args.scale = 2^-12*1.1;
    args.showImg = true;
    args.voltage = 140;

    for i=1:2:length(varargin)
        if(isfield(args,varargin{i}))
            eval(['args.' varargin{i} ' = varargin{i+1};']);
        else
            error(sprintf('unrecognized field for arguments in genBeamHarden: %s', varargin{i}));
        end
    end

    if(~isfield(args,'epsilon') || ~isfield(args,'iota'))
        [args.epsilon,args.iota]=readSpectrum('tungsten',args.voltage,0.05);
    end
end

