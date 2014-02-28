function [x, f, gF, Out]= flbfgsb(x, cT, AT, b, MT, lbd, ubd, nbd, opts)

%--------------------------------------------------------------------------
%
% L-BFGS method for solving
% min c'*x + 0.5* ||AT*x - b||^2, subject to l <= x <= u
%
% Copyright (c) 2008.  Zaiwen Wen, Wotao Yin
% Based on lbfgs.m by Liam Stewart
%--------------------------------------------------------------------------

% calculate objective value using the given mu_
    function val = calc_f()
        % nrm1x = norm(x,1);
        if isempty(MT); qx = 0.5*norm(res)^2; else qx = 0.5*res'*(MT*res); end
        val = cT*x + qx;  Out.nfe = Out.nfe+1;  %  it's safe to add an extra "sum" if mu a scalar,
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
        Out.nge = Out.nge+1; gF.Ax = Ax; gF.qx = qx; 
    end
% state var update: this function update status variables whenenve a new x is generated
    function x_update()
        % residual: used by calc_f and calc_g
        if isstruct(AT); xtmp=zeros(AT.n,1); xtmp(AT.col)=x; Ax = AT.A*xtmp; else Ax = AT*x; end
        res = Ax-b; 
    end
%     function x_update_ls()
%         % residual: used by calc_f and calc_g
%         Ax = Axp + stp*Ad; res=Ax-b;
%         nrmx = norm(x);
%     end
%     function x_prev_save()
%         xp = x; gp = g;  Axp = Ax; resp = res; fp = f; nrmxp = nrmx; 
%     end

n = length(x(:)); f = 0; g = zeros(n, 1);

res = zeros(length(b),1); fp = 0; qx = 0; xp = zeros(n,1); gp = zeros(n,1);  Ax = []; Axp = []; gF = []; 
Out.nfe = 0;  Out.nge = 0;

%--------------------------------------------------------------------------
% 
% if isempty(ubd) && isempty(lbd)
%     ubd = zeros(n,1); lbd = zeros(n,1); nbd = zeros(n,1);
% end
% if isempty(lbd)
%     lbd = zeros(n,1);
% end
% if isempty(nbd)
%     nbd = zeros(n,1);
% else
%     %          nbd(i) = 0 if x(i) is unbounded
%     %                   1 if x(i) has a lower bound (in lb(i))
%     %                   2 if x(i) has both lower (lb(i)) and upper (ub(i)) bounds
%     %                   3 if x(i) has an upper bound (ub(i))
%     idx = (
% end

assert(length(nbd) == n);
assert(length(ubd) == n);
assert(length(lbd) == n);

id = lbfgs_mex('init', n, x, lbd, ubd, nbd, opts);

maxits = opts.maxits;
Out.exitflag = 0;
it = 1; fp = 0; f = 0; fp_f_max = 0;
while 1
    [x,cmd] = lbfgs_mex('step', id, f, g);
    
    if cmd == 0                         % evaluate f, g
        x_update();    f = calc_f();  g = calc_g();
    elseif cmd == 1                     % iteration complete

        fp_f = abs(fp-f);
        %fprintf('it: %3d, f: %3.2e, fp: %3.2e, fp-f: %3.2e, fp_f_max: %3.2e, ||g||: %3.2e\n', it, f, fp, fp-f, fp_f_max, norm(g)); 
        if it >= maxits || abs( fp_f/max(abs([f,fp,1]))) <= min(opts.factr,1e-16) %|| fp_f <= eta_ffp*fp_f_max;
            Out.exitflag = 0;
            lbfgs_mex('stop', id);
            break;
        else
            it = it + 1;
            fp_f_max = max(fp_f_max, fp_f); fp = f;
        end
    elseif cmd == 2                     % converged
        Out.exitflag = 0;
        break;
    elseif cmd == 3                     % abnormal termination
        Out.exitflag = 1;
        break;
    elseif cmd == 4                     % error
        Out.exitflag = 2;
        break;
    else                                % unknown
        Out.exitflag = 3;
        break;
    end
end

lbfgs_mex('destroy', id);

end
