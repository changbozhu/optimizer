function [ out_y ] = dphi( obj_grad, xk, pk, alpha)
% Elvaluate gradient of phi(alpha) at the point alpha:
%           dphi(alpha) = nabla f(xk + alpha*pk)'pk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    AUTHER: Changbo Zhu
%    E-mail: changbozhu@outlook.com           
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xk=xk+alpha*pk;
    fgrad_xk = obj_grad(xk);
    out_y =dot(fgrad_xk,pk);
    return;

end

