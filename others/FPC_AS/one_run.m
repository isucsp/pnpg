function one_run
% this file runs a quick test on FPC_AS for solving the problem
%
%    min_x mu*|x|_1 + 0.5*|Ax-b|_2^2,
%
% where A is a m-by-n, and b=A*xs, where xs is a k-sparse vector. Because b
% contains no noise, a tiny mu of 1E-10 is used. This extremely small value
% of mu makes this problem numerically challegning.

% clc; 
clear all; clear classes; close all; 

mu=1e-10; 
seed = 200;                 % seed for random number generators in getData.m
n = 2^9;
delta = 0.5;                % m/n, m = round(delta*n)
rho = 0.3;                  % k/m, k = round(rho*m)
Ameth = 0;                  % see getData.m for help
xmeth = 1;                  %  "      "      "    "
% set noise level (mu depends on sigma1 and sigma2)
sigma1 = 0;    %- standard deviation of signal noise (added to xs)
sigma2 = 0;    %- standard deviation of meas. noise (added to b)
% problem size
m = round(delta*n); k = round(rho*m);
% initialization, get problem
[A,b,xs,xsn] = getData(m,n,k,Ameth,xmeth,sigma1,sigma2,seed); M = []; opts.xs = xs;
% call FPC_AS

opts.gtol = 1e-8;   % a termination option of FPC_AS; see manual
[x, Out] = FPC_AS(n,A,b,mu,M,opts); % call FPC_AS

% plot solutions
fig = figure(1); idn = 1:n;

epsIx = 0.1*min(abs(xs(xs~=0)));
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
title('recovered solution');
hold on;
xlabel('index'); ylabel('absolute value of entry, shown log-sccale');