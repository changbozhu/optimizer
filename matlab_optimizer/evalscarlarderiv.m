function [deriv, min_err] = evalscarlarderiv(phi, varargin)
% Evaluates derivative of the function phi at point 0 according to Ridders' method:
%
% Ridders, CJF. (1982). Accurate computation of F'(x) and F'(x) F''(x). Advances in Engineering
%     Software, 4(2), 75-6.
%
% INPUT:
%    phi           Function handle of a scalar real function of one real variable
%
% OUTPUT:
%    deriv         First derivative of phi at 0
%    min_err       Error estimate
%
% OPTIONS:
%    Optionally, the third argument of the function can be a structure containing further
%    settings for Ridder's method.
%
%    varargin{1}.init_eps       Initial  finite difference (default: 1)
%    varargin{1}.ratio          which is used to reduce eps on each step (default: 0.8)
%    varargin{1}.min_steps      Minimum number of steps in eps (default: 5)
%    varargin{1}.max_steps      Maximum number of steps in eps (default: 100)
%    varargin{1}.threshold      Terminate if last step worse than preceding
%                               by a factor of threshold
%                               (default: 1.5)
%    varargin{1}.tolerance      terminal computation if the fact error is
%                               less than tolerance (default: 0.001)
%
% --------------------------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    AUTHER: Changbo Zhu
%    E-mail: changbozhu@outlook.com           
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    warning off;
    % result
    deriv           = NaN;
    if nargin > 1
        options = varargin{1};
        options = evalderiv_options(options);
    else
        options = evalderiv_options();
    end
        
    % upacking parameters    
    min_err   = options.min_err;
    tolerance = options.tolerance;
    if isfield(options,'init_eps')
        init_eps = options.init_eps;
    end
        
    if isfield(options,'ratio')
        ratio = options.ratio;
    end
        
    if isfield(options,'min_steps')
        min_steps = options.min_steps;
    end
        
    if isfield(options,'max_steps')
        max_steps = options.max_steps;
    end
        
    if isfield(options,'threshold')
        threshold = options.threshold;
    end
    
    
    % Initialize recording matrix A by max_steps times max_steps
    A = NaN(max_steps);
    
    % Initialize finite difference step epsilon (eps)
    eps = init_eps;
    M = max_steps - 1;
    N = max_steps - 1;
    % Compute central difference formula at initial step
     
    A(1,1) = ( phi(eps) - phi(-eps) )/(2*eps);
    if isnan(A(1,1))
        disp("evalscarlarderiv.m-Warning: Richardson extrapolation A(1,1) is NaN")
    end
    % Loop to reduce eps 
    for m = 1:M
        n=0;
        % Generate a new step size eps_m = ratio^m * init_eps
        eps = ratio*eps;
        
        % Compute 2nd order approximation by central difference formula at 
        % this step 
        A(m+1,n+1) = ( phi(eps) - phi(-eps) )/(2*eps);
        % Use square of ratio for extrapolation because errors increase
        % quadratically with eps (here, of course, they decrease quadratically
        % because we're reducing eps ...)
        ratio_sq = ratio^2;
        coeff = ratio_sq;
        % Fill the current row using Richardson extrapolation
        for n = 1:m 
            % Richardson extrapolation
            A(m+1,n+1) = ( A(m+1,n) - coeff*A(m,n) )/(1 - coeff);
            if isnan(A(m+1,n+1))
                disp(['evalscarlarderiv.m-Warning: Richardson extrapolation A(' ...
                    num2str(m+1) ',' num2str(n) ') is NaN'])
            end
            % Increment extrapolation factor
            coeff = coeff*ratio_sq;
            
            % Error on this trial is defined as the maximum absolute difference
            % to the extrapolation parents
            err_mn = max(abs(A(m+1,n+1)-A(m+1,n)),abs(A(m+1,n+1) - A(m,n)));
            if err_mn < min_err
                min_err = err_mn;
                deriv = A(m+1,n+1);
            end
        end

        % Stop if errors start increasing (to be expected for very small
        % values of eps
        if m > (min_steps - 1) && abs( A(m+1,n+1) - A(m,n) ) > (threshold*min_err)...
           && min_err < tolerance
            return
        end
    end
end