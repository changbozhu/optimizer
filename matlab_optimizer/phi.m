function [ out_y ] = phi(obj_func, xk, pk, alpha )
% Definition of optimised objective function phi:
%           phi(alpha) = f( xk+alpha*pk )
% where f is the optimized objective function obj_func
    xk = xk + alpha * pk;
    out_y=obj_func(xk);
    return;
end

