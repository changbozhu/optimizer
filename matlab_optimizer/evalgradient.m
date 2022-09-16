function [grad, min_err] = evalgradient(f, x, varargin)
% Calculates the gradient of the function f at point x according to Ridders' method:
%
% Ridders, CJF. (1982). Accurate computation of F'(x) and F'(x) F''(x). Advances in Engineering
%     Software, 4(2), 75-6.
%
% INPUT:
%    f             Function handle of a real function of n real variables which are passed as
%                  *one* vector with n elements
%    x             Point at which to differentiate f
%
% OUTPUT:
%    gradf         Gradient of f at x (a column vector)
%    min_err           Error estimates (a column vector)
%
% OPTIONS:
%    Optionally, the third argument of the function can be a structure containing further
%    settings for Ridder's method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    AUTHER: Changbo Zhu
%    E-mail: changbozhu@outlook.com           
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
    warning off;
    if any(isnan(x)) || any(isinf(x))
        error('x is not a numerical vector!')
    end
    d = length(x);
    x = reshape(x,d,1);
    I = eye(d);
    grad   =  NaN(d,1);
    min_err =  realmax * ones(d,1);
    % Biuld struct of configuration for evaluating gradient of function f.    
    if nargin > 2
        options = varargin{1};
        options = evalgradient_options(options);
    else
        options = evalgradient_options();
    end    
    
    % Check if f and x match
    try
        FrE = f(x);
        if isinf(FrE) || isnan(FrE)
            error(['evalgradient.m-FrE=' num2str(FrE)])
        end
    catch err
        error('Function cannot be evaluated at differentiation point');
    end
    
    if options.mpi == true
        % Loop through each component of variable x by using parfor     
        %startmatlabpool;
        parfor i = 1:d
            ei = I(:,i);
            % Construct filehandle to be passed to eval
            phi_eps = @(eps) phi(eps,f,x,ei);        
            % Calculate derivative
            [grad(i), min_err(i)] = evalscarlarderiv(phi_eps,options);            
        end
        %closematlabpool;
    else
        for i = 1:d
            ei = I(:,i);
            % Construct filehandle to be passed to eval
            phi_eps = @(eps) phi(eps,f,x,ei);
            % Calculate derivative
            [grad(i), min_err(i)] = evalscarlarderiv(phi_eps,options);
        end
    end
    
    
end

function y = phi(eps,f,x,p)
    d = length(x);
    x = reshape(x,d,1);
    p = reshape(p,d,1);
    y = f(x + eps*p);
end
