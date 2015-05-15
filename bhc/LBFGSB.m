classdef LBFGSB < handle
    properties
        Ie;
        func;
        thresh=1e-5;
        cost;
        maxItr=1e3;
        opts

        zmf;
        warned=false;
    end
    methods
        function obj = LBFGSB(Ie,maxItr)
            obj.Ie=Ie;
            obj.opts.maxIts=maxItr;
            obj.opts.maxTotalIts=5e5;
            obj.opts.printEvery=inf;
        end
        function [Ie,numItr]=main(obj)
            obj.opts.maxIts=obj.maxItr;
            obj.opts.x0=obj.Ie;
            obj.opts.factr=obj.thresh/eps;
            [Ie,cost,info]=lbfgsb(obj.func,zeros(size(obj.Ie)),inf*ones(size(obj.Ie)),obj.opts);
            obj.Ie=Ie;
            obj.cost=cost;
            numItr=info.iterations;
        end
    end
end
