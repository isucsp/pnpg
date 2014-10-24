function [x, f, gF, Out]= fminLBFGS(x, cT, AT, b, MT, opts, varargin)

%--------------------------------------------------------------------------
%
% L-BFGS method for solving
% min c'*x + 0.5* ||AT*x - b||^2
%
% Copyright (c) 2008.  Zaiwen Wen, Wotao Yin
%
%--------------------------------------------------------------------------

%     function y = ProdA(x)
%        y = A*x;  if smax ~= 1; y=y/smax; end; Out.nfe = Out.nfe+1;
%     end
%     function y = ProdAT(x)
%        if isa(A,'A_operator'), y = AT_times_x(A,x); else y = A'*x; end;
%        if smax ~= 1; y=y/smax; end; Out.nge = Out.nge+1;
%     end

% calculate objective value using the given mu_
    function val = calc_f()
        if isempty(MT); qx = 0.5*norm(res)^2; else qx = 0.5*res'*(MT*res); end;  val = cT*x + qx;  
    end
% calculate gradient of the fidelity term
    function grad = calc_g()
        if isstruct(AT);
            if isempty(MT); gF.g = (AT.A')*res;  else gF.g = (AT.A')*(MT*res);  end
            grad = cT' + gF.g(AT.col); gF.res = res;
        else
            if isempty(MT); gF.gT = res'*AT;   else gF.gT = res'*MT*AT;   end
            grad = (cT + gF.gT)'; gF.res = res;
        end
        nrmG = norm(grad); Out.nge = Out.nge+1; gF.Ax = Ax; gF.qx = qx; 
    end
% state var update: this function update status variables whenenve a new x is generated
    function x_update()
        if isstruct(AT); xtmp=zeros(AT.n,1); xtmp(AT.col)=x; Ax = AT.A*xtmp; else Ax = AT*x; end
        res = Ax-b; nrmx = norm(x); Out.nfe = Out.nfe+1;  
    end
% residual: used by calc_f and calc_g
    function x_update_ls()
        Ax = Axp + stp*Ad; res=Ax-b;  nrmx = norm(x);
    end
    function x_prev_save()
        xp = x; gp = g;  Axp = Ax; resp = res; fp = f; nrmxp = nrmx; gFp = gF; 
    end

n = length(x);
opts.memo =[];
% termination rule
opts = AddLowerBoundedNumberOption(opts, 'debug', 0, 0, 1,'print debug information');
opts = AddLowerBoundedNumberOption(opts, 'maxitr', 1000, 1, 1e5, 'max number of iterations');
opts = AddLowerBoundedNumberOption(opts, 'xtol', 1e-8, 0, 1, 'Tolerance on norm(x - xp)/norm(xp)');
opts = AddLowerBoundedNumberOption(opts, 'gtol', 1e-8, 0, 1, 'Tolerance on norm of gradient');

opts = AddLowerBoundedNumberOption(opts, 'm', 5, 1, 100, 'storage number of L-BFGS');
opts = AddLowerBoundedNumberOption(opts, 'rho1', 1e-4, 0, 1, 'parameters for control the linear approximation in line search');
opts = AddLowerBoundedNumberOption(opts, 'rho2', 0.9, 0, 1, 'parameters for control the linear approximation in line search');
parsls.ftol = opts.rho1;
parsls.gtol = opts.rho2;

%-------------------------------------------------------------------------------
% copy parameters
xtol = opts.xtol;
gtol = opts.gtol;
m = opts.m;

res = zeros(length(b),1); fp = 0; qx = 0; xp = zeros(n,1); gp = zeros(n,1);  Ax = []; Axp = []; gF = []; % resp = zeros(m,1);
Out.nfe = 0;  Out.nge = 0;
istore = 0; pos = 0;  ppos = 0;  perm = [];

%% Initial function value and gradient
% prepare for iterations
x_update(); f=calc_f(); g=calc_g();

% set up storage for L-BFGS
% if ~isfield(opts, 'istore'); istore = 0; else istore = opts.istore; end

SK = zeros(n,m);		  % S stores the last ml changes in x
YK = zeros(n,m);		  % Y stores the last ml changes in gradient of f.
STY = zeros(m);
YTY = zeros(m);

% if istore > 0
%     SK = opts.SK(:,1:istore);		  % S stores the last ml changes in x
%     YK = opts.YK(:,1:istore);		  % Y stores the last ml changes in gradient of f.
%     STY = opts.STY(1:istore,1:istore);
%     YTY = opts.YTY(1:istore,1:istore);
% end


% Print iteration header if debug == 1
if (opts.debug == 1)
    fid = 1;
    fprintf(fid, '\n----------- BFGS Method with CSRCH Line Search ----------- \n');
    fprintf(fid, '%4s \t %10s \t %10s \t %10s \t %10s\n', 'Iter', 'dt', 'f', 'CritdiffX', '||G||' );
    fprintf(fid, '%4d \t %4.3e \t %4.3e \t %4.3e \t %4.3e\n', 0, 0, f, 0, nrmG);
    %OUTPUT
end

% main loop
for iter = 1:opts.maxitr

    % store old point
    x_prev_save(); 

    % compute search direction
    % if the first iteration
    if istore == 0
        d = -g; stp = 1;
    else
        gamma = STY(ppos, ppos)/YTY(ppos, ppos);
        Rinv = inv(STY(1:ppos,1:ppos));
        d = -( gamma*g + [SK(:,perm2), gamma*YK(:,perm2)] * ...
             ([ Rinv'*( diag(diag(STY(1:ppos,1:ppos))) + gamma* YTY(1:ppos,1:ppos)) * Rinv , -Rinv';
              - Rinv, zeros(ppos)]*[SK(:,perm2)'*g; gamma*(YK(:,perm2)'*g)]));
        stp = 1;
    end

    x = x + d;     x_update();    f = calc_f(); % new g is not needed at this moment
    % do exact line search since the objective function is quadratic
    Ad = Ax - Axp; c1 = Ad'*Ad; c2 = Ad'*resp + cT*d;
    stp = -c2/c1; if stp < 0 || isnan(stp) || isinf(stp); stp = 0.1; end
    x = xp + stp*d; f = 0.5*stp^2*c1 + stp*c2 + fp; 
    
    if f < 0;
       x = xp; f = fp; gF = gFp;
       Out.msg = 'negative function value';
       Out.iter = iter;
       Out.nge = Out.nfe;
       return;
    end
    
    x_update_ls(); g=calc_g(); 
%     if stp > 0 && ~isnan(stp) && ~isinf(stp); 
%         x = xp + stp*d; f = 0.5*stp^2*c1 + stp*c2 + fp; 
%     else
%         if (opts.debug == 1); fprintf('Exact stp is: %4.3e, either less than zero, or inf or nan \n',stp); end
%         % begin line search, must set "work.task = 1"
%         workls.task =1;     deriv = c1 + c2;
%         % call line search, reverse communication
%         while 1
%             [stp, f, deriv, parsls, workls] = ls_csrch(stp, f, deriv , parsls , workls);
%             % Evaluate the function and the gradient at stp
%             if (workls.task == 2)
%                 % calculate g, f,
%                 x = xp + stp*d;
%                 f = 0.5*stp^2*c1 + stp*c2 + fp;
%                 deriv = stp*c1 + c2;
%             else  % exit
%                 break
%             end
%         end
%         if (opts.debug == 1); fprintf('Inexact stp is: %4.3e\n',stp); end
%     end
%     x_update_ls(); g=calc_g();
    
    % s = x - xp = stp*d;  ==> ||s|| = stp*||d||
    nrms = stp*norm(d);
    % compute stopping 
    CritdiffX = nrms/max(nrmxp,1);
    
    % now, update normG
    Out.nrmg =  nrmG;
    
    if (opts.debug == 1)
       fprintf(fid, '%4d \t %4.3e \t %4.3e \t %4.3e \t %4.3e\n', iter, stp, f, CritdiffX, nrmG);
    end
    
    fp_f = abs(fp-f);
    %if (CritdiffX < xtol) || (nrmG < gtol*max(1,nrmx)) || fp_f/max(abs([f,fp,1])) <= xtol*1e-10
    if (nrmG < gtol*max(1,nrmx)) || fp_f/max(abs([f,fp,1])) <= xtol*1e-10
        
        Out.msg = 'converge';
        Out.iter = iter;
        Out.nge = Out.nfe;
        return;
    end
    
    %----------------------------------------------------------------------
    % save for L-BFGS
    ygk = g-gp;		%save change in g for L-BFGS.
    s = x-xp;

    %Check to save s and y for L-BFGS.
    if ygk'*ygk>1e-20
        istore = istore + 1;
        pos = mod(istore, m); if pos == 0; pos = m; end;    ppos = pos; 
        YK(:,pos) = ygk;  SK(:,pos) = s;

        if istore <= m
            perm = [perm, pos]; perm2 = perm;
            % STY is S'Y, upper triangular
            STY(1:istore,istore)  = SK(:,1:istore)'*ygk;
            % YTY is Y'Y
            YTY(1:istore,istore)  = YK(:,1:istore)'*ygk;
            YTY(istore,1:istore)  = YTY(1:istore,istore)';
        else
            % update YK, SK
            ppos = m; perm = [perm(m), perm(1:m-1)]; perm2 = [perm2(2:m), perm2(1)];
            % update STY, YTY, first move the old part to upper
            STY(1:end-1, 1:end-1) = STY(2:end, 2:end);
            YTY(1:end-1, 1:end-1) = YTY(2:end, 2:end);

            % then update the last column or row
            STY(perm,end)  = SK(:,:)'*ygk;
            YTY(perm,end)  = YK(:,:)'*ygk;
            YTY(end,:)  = YTY(:,end)';
        end
    end    
    
end

Out.msg = 'Exceed max iteration';
Out.iter = iter-1;
end