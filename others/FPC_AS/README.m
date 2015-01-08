% FPC_AS: A MATLAB Toolbox for solving l1 minimization
%
% -------------------------------------------------------------------------
% 1. Introduction
%
% Thank you for downloading the FPC_AS toolbox! FPC_AS is a package for
% solving
%   min mu*||x||_1 + 0.5*||Ax - b||_M^2,
% using an active set continuation algorithm using shrinkage and 
% sub-optimization.
%
% A detailed explanation of the toolbox is given in FPC_Active_Manual.pdf.
% 
% -------------------------------------------------------------------------
% 2. Optional external package
%
% FPC_AS is prepackaged with all required dependencies, but if you wish to
% use the code L-BFGS-B to solve sub-optimization with bounded constraints,
% then you need to install:
% 
%     MATLAB wrapper for L-BFGS-B   
%       http://www.cs.toronto.edu/~liam/software.shtml
%   
%     L-BFGS-B 
%       http://www.ece.northwestern.edu/~nocedal/lbfgsb.html
%           
% -------------------------------------------------------------------------
% 3. Installation and Setup
% 
%   3.1) Remove any old version of FPC_AS
%   3.2) unzip FPC_AS.zip. This should create the structure
%       /FPC_AS, /FPC_AS/src, 
%       /FPC_AS/prob_gen, /FPC_AS/prob_gen/classes 
%   3.3) Add these folders to your MATLAB path
%   3.4) Add the path of the MATLAB wrapper for L-BFGS-B (if installed)
%       to your MATLAB path
%
% -------------------------------------------------------------------------
% 4. Quick Start:
% 
% [x, Out] = FPC_AS(n,A,b,mu);
% 
% n the dimension of x, the nubmer of columns of A
%
% A either an explicit m x n matrix or an A_operator object representing
%   a matrix implicitly. When A is relatively big and/or operations A*x 
%   and A.'*x can be computed much faster through certain means, it is 
%   recommend that A is created as an A_operator object, the implementation
%   of which is provided with the solver. 
%
%   To create an A_operator object, two functions or function handles for
%   computing A*x and A.'*x, respectively, must be given. Suppose they
%   are AA and AT,
%       AA: a function handle such that AA(x) = A*x,
%       AT: a function handle such that AT(x) = A.*x.
%   Then A can be created as an A_operator by
%       A = A_operator(@(x) AA(x), @(x) AT(x));
%
%   An example for A being an implicit matrix that performs a discrete 
%   cosine transform (DCT) and returns the subset of the DCT result 
%   corresponding to omega is
%
%     function y=pdct(x,picks); y=dct(x); y=y(picks); end
%     function y=pidct(x,n,picks); y=zeros(n,1); y(picks)=x; y=idct(y); end
%     A = A_operator(@(x) pdct(x,omega), @(x) pidct(x,n,omega)); 
%
% b must be an m x 1 vector, the row rank of A must equal to m
%
% For more detailed usage, please see the manual FPC_Active_Manual.pdf or
% the example driver file one_run.m
%
% -------------------------------------------------------------------------
% 5. The Test Problems
% A collection of test problems is stored in "prob_gen"
%
%
% -------------------------------------------------------------------------
% 6. License
% 
% See License.txt for the license of this program
%
% -------------------------------------------------------------------------
% 7. The Authors
%
% We hope that FPC_AS is useful for your application.  If you have
% any bug reports or comments, please feel free to email one of the
% toolbox authors:
%
%   Wotoa Yin,  wotao.yin@rice.edu
%   Zaiwen Wen, zw2109@columbia.edu
%
% Enjoy!
% Wotao and Zaiwen
%
%

