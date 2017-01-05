function [D,iter,fun_all]=denoise_bound_mod(Xobs,lambda,l,u,pars)
% optimized by renliang gu
%
%This function implements the FISTA method for TV-based denoising problems
%
% Based on the paper
% Amir Beck and Marc Teboulle, "Fast Gradient-Based Algorithms for Constrained
% Total Variation Image Denoising and Deblurring Problems"
% -----------------------------------------------------------------------
% Copyright (2008): Amir Beck and Marc Teboulle
% 
% FISTA is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
% INPUT
% Xobs ..............................an observed noisy image.
% lambda ........................ parameter
% l ..................................... lower bound on the pixels' values
% u ..................................... upper bound on the pixels' values
% pars.................................parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.epsilon ..................... tolerance for relative error used in
%                                                       the stopping criteria (Default=1e-4)
% pars.print ..........................  1 if a report on the iterations is
%                                                       given, 0 if the  report is silenced
% pars.tv .................................. type of total variation
%                                                      penalty.  'iso' for isotropic (default)
%                                                      and 'l1' for nonisotropic
%  
% OUTPUT
% D ........................... The solution of the problem 
%                                            min{||X-Xobs||^2+2*lambda*TV(X
%                                            ) : l <= X_{i,j} <=u} 
% iter .............................  Number of iterations required to get
%                                            an optimal solution (up to a tolerance)
% fun_all ......................   An array containing all the function
%                                             values obtained during the
%                                             iterations

%Define the Projection onto the box
if((l==-Inf)&&(u==Inf))
    project=@(x)x;
elseif (isfinite(l)&&(u==Inf))
    project=@(x)max(x,l);
elseif (isfinite(u)&&(l==-Inf))
    project=@(x)min(x,u);
elseif ((isfinite(u)&&isfinite(l))&&(l<u))
    project=@(x) max(min(x,u),l);
else
    error('lower and upper bound l,u should satisfy l<u');
end

% Assigning parameres according to pars and/or default values
flag=exist('pars','var');
if (flag&&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
if (flag&&isfield(pars,'epsilon'))
    epsilon=pars.epsilon;
else
    epsilon=1e-4;
end
if(flag&&isfield(pars,'print'))
    prnt=pars.print;
else
    prnt=1;
end
if(flag&&isfield(pars,'tv'))
    tv=pars.tv;
else
    tv='iso';
end

[m,n]=size(Xobs);
D=zeros(m,n);
if(isfield(pars,'P1') && isfield(pars,'P2'))
    P1=pars.P1; P2=pars.P2;
else
    [P1,P2]=Ltrans(D);
end

switch tv
    case 'iso'
        A=[P1;zeros(1,n)].^2+[P2,zeros(m,1)].^2;
        A=sqrt(max(A,1));
        P1=P1./A(1:m-1,:);
        P2=P2./A(:,1:n-1);
    case {'l1', '1d'}
        P1=max(min(P1,1),-1);
        P2=max(min(P2,1),-1);
    otherwise
        error('unknown type of total variation. should be iso or l1');
end
R1=P1; R2=P2;

% debug=true;
f=@(x,u) 0.5*norm(x-Xobs,'fro')^2+u*tlv(x,tv);

tk=1;
tkp1=1;
count=0;
i=0;

fval=inf;
fun_all=[];
if(prnt)
    fprintf('***********************************\n');
    fprintf('*Solving with FGP/FISTA**\n');
    fprintf('***********************************\n');
    fprintf('#iteration  function-value  relative-difference\n');
    fprintf('---------------------------------------------------------------------------------------\n');
end
while((i<MAXITER)&&(count<3))
    fold=fval;
    %%%%%%%%%
    % updating the iteration counter
    i=i+1;
    %%%%%%%%%
    % Storing the old value of the current solution
    Dold=D;
    %%%%%%%%%%
    %Computing the gradient of the objective function
    Pold1=P1;
    Pold2=P2;
    tk=tkp1;
    C=Xobs-lambda*Lforward(R1,R2);
    D=project(C);
    [Q1,Q2]=Ltrans(D);

    %%%%%%%%%%
    % Taking a step towards minus of the gradient
    P1=R1+Q1/(8*lambda);
    P2=R2+Q2/(8*lambda);

    %%%%%%%%%%
    % Peforming the projection step
    switch tv
        case 'iso'
            A=[P1;zeros(1,n)].^2+[P2,zeros(m,1)].^2;
            A=sqrt(max(A,1));
            P1=P1./A(1:m-1,:);
            P2=P2./A(:,1:n-1);
        case {'l1','1d'}
            P1=max(min(P1,1),-1);
            P2=max(min(P2,1),-1);
            %     P1=P1./(max(abs(P1),1));
            %     P2=P2./(max(abs(P2),1));
        otherwise
            error('unknown type of total variation. should be iso or l1');
    end

    %%%%%%%%%%
    %Updating R and t
    tkp1=(1+sqrt(1+4*tk^2))/2;

    R1=P1+(tk-1)*(P1-Pold1)/tkp1;
    R2=P2+(tk-1)*(P2-Pold2)/tkp1;

    re=norm(D-Dold,'fro')/norm(D,'fro');
    if (re<epsilon)
        count=count+1;
    else
        count=0;
    end

    fval=(-norm(C-D,'fro')^2+norm(C,'fro')^2)/2;
    if (fval>fold)
        tkp1=1;
    end

    if(nargout==3)
        fun_all=[fun_all;fval];
        if(prnt)
            fprintf('%7d %10.10g %10.10g  %19g',i,fval,norm(D-Dold,'fro')/norm(D,'fro'),f(D,lambda));
            if (fval>fold) fprintf('  *\n'); else fprintf('   \n'); end
        end
    end
end
iter=i;
end

