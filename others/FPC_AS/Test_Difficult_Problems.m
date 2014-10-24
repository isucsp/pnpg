function Test_Difficult_Problems()



% clc;
clear classes; close all; 
% addpath(genpath('./'));

mu=1e-10;


%% options for FPC_AS; see manual
% opts.mxitr = 500;
opts.record = true;
opts.record = 0;

% The default of opts.zero has been changed to 1e-10. 
% opts.zero = 1e-8;  
% opts.gtol = 1e-14;
opts.gtol = 1e-12;
opts.gtol_scale_x = 1e-16;
% opts.PrintOptions = 0;
% opts.scale_A = 0;
% opts.dynamic_zero = 0;
%
% opts.ls_meth = 'nmls';
opts.ls_meth = 'Hybridls';
% opts.sub_opt_meth = 'lbfgsb';
opts.sub_opt_meth = 'lbfgs';
opts.sub_opt_meth = 'pcg';

Problist = {'CaltechTest1', 'CaltechTest2', 'CaltechTest3', 'CaltechTest4', ...
    'Ameth6Xmeth2n1024m512k150seed200', 'Ameth6Xmeth2n1024m512k151seed200',...
    'Ameth6Xmeth2n1024m512k152seed200', 'Ameth6Xmeth2n1024m512k153seed200', ...
    'Ameth6Xmeth2n1024m512k154seed200', 'Ameth6Xmeth6n1024m512k154seed200'};


% Problist = {'Ameth6Xmeth2n1024m512k150seed200'};
% Problist = {'Ameth6Xmeth2n1024m512k151seed200'};
% % Problist = {'Ameth6Xmeth2n1024m512k152seed200'};
% % Problist = {'Ameth6Xmeth2n1024m512k153seed200'};
% Problist = {'Ameth6Xmeth2n1024m512k154seed200'};
% Problist = {'Ameth6Xmeth6n1024m512k154seed200'};

% Problist = {'CaltechTest1', 'CaltechTest2', 'CaltechTest3', 'CaltechTest4'};
% Problist = {'CaltechTest1'};
% Problist = {'CaltechTest2'};
% Problist = {'CaltechTest3'};
% Problist = {'CaltechTest4'};

% Test problems from Candes
for di = 1:length(Problist)
    fprintf('\n\nTest problem: %s', Problist{di});
    % load problem
    switch Problist{di}
        case { 'BadEx1'}
            load('Data_p=1_k=31_rep=3.mat','m','n','Ameth','A','b','xs');
            alpha = 0.5; full = true; sig1=0; sig2=0; opts.xs = xs;
            [M,mu,A,b,sig,kap,tau,M12] = getM_mu(full,mu,m,n,Ameth,A,b,sig1,sig2,alpha);
         case  {'CaltechTest1', 'CaltechTest2', 'CaltechTest3', 'CaltechTest4', ...
                 'Ameth6Xmeth2n1024m512k154seed200','Ameth6Xmeth6n1024m512k154seed200', ...
                 'Ameth6Xmeth2n1024m512k150seed200', 'Ameth6Xmeth2n1024m512k151seed200', ...
                 'Ameth6Xmeth2n1024m512k152seed200', 'Ameth6Xmeth2n1024m512k153seed200'}            
            
            clear b Omega n xs x A 
            load( Problist{di}, 'b','Omega','n','xs');
            %set up problem
            m = length(b); A = A_operator(@(x) pdct(x,1,n,Omega), @(x) pdct(x,2,n,Omega)); M = []; opts.xs = xs;
%             m = length(b); A = dctmtx(n); A = A(Omega,:); M = []; opts.xs = xs;
            
    end


    epsIx = 0.1*min(abs(xs(xs~=0)));

    %---------------------------------------------------------------
    [x, Out] = FPC_AS(n,A,b,mu,M,opts);
    %---------------------------------------------------------------
    % plot

    fig = figure(di); idn = 1:n;

    if isempty(xs);
        semilogy(idn, abs(x),'.');  legend('x',0);
    else
        idxs = (xs~=0);
        semilogy(idn(idxs), abs(x(idxs)),'*', idn(~idxs), abs(x(~idxs)),'.', 'MarkerSize',6);
        hold on; semilogy(idn, abs(xs), 'dr', 'MarkerSize',6);
        legend('x on T', 'x on T^c','x^*');
        ypos = min(epsIx,1e-3); semilogy(idn, ypos*ones(n,1));
        text( n*1.01 , ypos*1.01, sprintf('Level: %2.1e',ypos), 'HorizontalAlignment','left')
        hold off;
    end;

    axis([0,n*1.25, 0, max(max(abs(x)))*2]);
    title([Problist{di} ' :recovered solution']);
    hold on;
    xlabel('index'); ylabel('absolute value of entry, shown log-sccale');

end


end

%--------------------------------------------------------------------------
% Note 1:
% The nummber of A'*x in which A is an explicit matrix might be different
% from the number of A'*x in which A is an A_operator since in the later
% case the gradient from sub-optimization is re-used in the stage of
% shrinkage.
%
%
%
% Note 2:
% The numerical results reported in the paper used the thresholding value
% "0.1*min(abs(xs(xs~=0)))" to compute the pair (sgn, miss, over) which
% might be different from the pair (sgn, miss, over) reported by FPC_AS 
% since default thresholding value in FPC_AS.m is opts.zero = 1e-10 
%

