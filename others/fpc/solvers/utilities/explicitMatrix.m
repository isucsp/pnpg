function y = explicitMatrix(A,x,mode)

% Function that can be used to create anonymous function A(x,mode) for 
% explicit matrices.

switch mode
    case 1
        y = A*x;
    case 2
        y = A'*x;
    otherwise
        error('Unknown mode passed to explicitMatrix in fpc.m');
end

end % explicitMatrix