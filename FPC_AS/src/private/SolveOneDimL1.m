function  [stp, f, x, Out] = SolveOneDimL1(n, stp1, stp2, x, d, mu, c1, c2, pars)

%
% min phi(a) = 0.5*c1*a^2 + c2*a + sum(mu.*|x+a*d|_1), subjec to a in [stp1, stp2]
% Author:  Zaiwen Wen
% Date: July, 5, 2008

    function f = feval
        xn = x+stp*d; f = 0.5*c1*stp^2 + c2*stp + sum(mu.*abs(xn));
    end

%% Step 0. Test if x_new = 0 ?
eps = 1e-16;

% check if stp = 1 optimal or not
stp = stp2; f = feval();
mud = mu.*d; 

%--------------------------------------------------------------------------
% A comment:
% const = opts.c*(gp'*dx);
% %const = opts.c*(gp'*dx + sum(mu.*(abs(x) - abs(xp))));
% If SolveOneDimL1 is called after f(x+stp2*d) > C + alpha*const 
% and if stp2 is optimal, is this a contradiction? 
% This case really happened. I am confused.
%--------------------------------------------------------------------------

% let g denote the subgradient of phi(a) at stp2
% if stp2 is optimal, then g <= 0
% phi(a) has terms like |x+a*d| and x+a*d might be zero
% The subgradient of the nonzero parts should be g1 below
%  case 1: g1 is less than zero
%  case 2: g1 is less than sum(abs(mu*d(E))) which is corresponding to the zero component
E = (abs(xn) <= eps); 
g1 = stp*c1+c2 + sum(mud(~E).*sign(xn(~E))); g2 = sum(abs(mud(E)));
if  g1 - g2 < 1e-6
    stp = 1;     x = xn;     Out.exit = 1; Out.mesg = 'stp2 is optimal';
    if pars.debug == 1; fprintf('stp=%4.3e is optimal, exit\n', stp2); end
%     fprintf('\n\n stp2 is optimal, g1: %4.3e, g2: %4.3e exit\n\n', g1, g2); 
    return;
end

% stp = -c2/c1; if stp < stp1; stp = stp2; end

%% Step 1. Initialization
E = (abs(d) <= eps); d(E) = 0;
IP = (d > eps); IN = (d < -eps); I = IP | IN;
beta = zeros(n,1); beta(I) = -x(I)./d(I);


%% Step 2. Compute lower and upper bounds for stp

stp_list = [beta( beta > stp1 & beta < stp2); stp1; stp2];
stp_list = sort(stp_list);
kappa =  length(stp_list);


%% Step 4. Compute stp
if pars.debug == 1; fprintf('stp: [%4.3e ,\t  %4.3e], # of intervals: %d\n', stp1, stp2, kappa-1); end

% main loop
Out.exit = 0;

% the "di" th interval means: [stp_list(di-1),  stp_list(di)]
% the set of intervals has not been marked yet
% a logical vector, if ind_mark(di) = true, the "di"th has not been
% searched yet, otherwise, it has been
%
ind_mark = true(kappa,1);  ind_mark(1) = false;

di = kappa;
% find an interval which includes the initial estimation of stp
if  (stp >= stp2 ||  stp <= stp1)
    % we begin with the last interval
    di =  kappa;
    stp = stp_list(kappa);
else
    di =  find(ind_mark & (stp_list>=stp)  , 1, 'first');
    % if di is empty
    if isempty(di)
        if pars.debug == 1
            fprintf('cannot find in the above\n');
        end
        di =  find(ind_mark & (stp_list<=stp)  , 1, 'last');
    end
end

% if there is no stp in [stp1, stp2], return the best value
Out.itr = 0; f_best = f; stp_best = stp;
for indi = kappa-1:-1:1 % only need to search kappa-1 intervals
    % take a interval
    Out.itr = Out.itr + 1;
    lb = stp_list(di-1);     ub = stp_list(di);
    % delete the interval which has been marked
    ind_mark(di) = false;

    % if lb == ub, the length o this interval is zero, skip it
    if abs(ub-lb)<=eps
        % pick the last element from the list
        di = find(ind_mark, 1, 'last');
        if isempty(di) || di == 1
            error(['Error: the next di is not available or di =1,', ...
                'therefore there is no interval']);
        end
        continue;
    end
    
    alpha = x + lb*d;   beta = x + ub*d; 
    ind = (alpha >= eps & beta >= eps); rho = sum(mud(ind));
    ind = (alpha < eps & beta < eps); rho = rho - sum(mud(ind));

    stp = -(c2+rho)/c1;
    f = feval();
   
    if pars.debug == 1
        fprintf('di:%d \t lb:%4.3e \t ub:%4.3e \t stp:%4.3e \t f:%4.3e \t stpb:%4.3e \t fb:%4.3e\n', di, lb, ub, stp, f, stp_best, f_best);
    end
    
    if  (stp <= ub && stp >= lb) || abs(stp-ub)< eps ...
                || abs(stp-lb)< eps
        x = xn;         Out.exit = 1;   Out.mesg = 'optimal solution';
        if pars.debug == 1; fprintf('optimal stp=%4.3e\n', stp); end
        return;
    elseif f_best > f && stp >= stp1 && stp <= stp2
        f_best = f; stp_best = stp;
    end
    
    if Out.itr == 5
        stp = stp_best; f = f_best; x = x + stp*d; Out.exit = 10; Out.mesg = 'exceed max iterations'; 
        if pars.debug == 1; fprintf('exceed max iterations\n', stp); end
        return;
    end
    
    % now, determine the index of the next interval from the intervals have not
    % been marked so far. We want stp is enclosed in this interval. We
    % first try to find a nearest interval greater than stp. If failed, we try
    % to find a nearest interval smaller than stp

    %----------------------------------------------------------------------
    di =  find(ind_mark & (stp_list>=stp)  , 1, 'first');
    % if di is empty
    if isempty(di)
        di =  find(ind_mark & (stp_list<=stp)  , 1, 'last');
    end
    %----------------------------------------------------------------------

    if isempty(di) || di == 1
        %error(['Error: the next di is not available or di =1, therefore there is no interval']);
        break;
    end

end

stp = stp_best; f = f_best; x = x + stp*d; Out.exit = -1; Out.mesg = 'no optimal solution is found, return the best so far';
% fprintf(' No stp is found in [%4.3e, %4.3e], return the best value so far\n', stp1, stp2);

% error(' \n  No stp is found, Debug this ! \n\n');
%determine if a stp  is found or not
if Out.exit == 0 && pars.debug == 1
    fprintf(' No stp is found in [%4.3e, %4.3e], return the best value so far\n', stp1, stp2);
    fprintf(' \n Debug this ! \n\n'); 
%     error(' \n Debug this ! \n\n');
end


end
