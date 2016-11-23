classdef Proximal < handle
    % interface class that solves: 
    % minimize  0.5*||x-a||_2^2+u*r(x)
    % where a and u are given when calling denoise function
    properties
        regFunc;    % penalty function handle r(x)
        iterative=false;
        steps;          % number of iteration have run
        maxItr=100;
        thresh=1e-9;
    end
    properties (Access = protected)
        init
        proximalOp;
        x
    end
    properties (Dependent)
    end
    methods
        function obj = Proximal(op,regf)
            obj.init=0;
            if(exist('op','var'))
                obj.proximalOp=op;
                if(nargin(op)>2) % That means proximal operator is iterative
                    obj.iterative=true;
                end
            end
            % In the future, need to transfer the ownership from alphaStep to Proximal
            if(exist('regf','var'))
                obj.regFunc=regf;
            end
        end
        function setInit(obj)
            obj.init=obj.x;
        end
        function set.maxItr(obj,maxItr)
            obj.maxItr=maxItr;
        end
        function set.thresh(obj,thresh)
            obj.thresh=thresh;
        end
        function x=denoise(obj,a,u)
            switch(nargin(obj.proximalOp))
                case 1
                    x=obj.proximalOp(a);
                case 2
                    x=obj.proximalOp(a,u);
                otherwise
                    [x,obj.steps]=obj.proximalOp(a,u,obj.thresh,obj.maxItr,obj.init);
            end
            obj.x=x;
        end
        function y=getPenalty(obj,x)
            if(~exist('x','var'))
                y=obj.regFunc(obj.x);
            else
                y=obj.regFunc(x);
            end
        end
    end
    methods (Access = private)
    end
end

