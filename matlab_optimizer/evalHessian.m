 function [Hessian, min_err] = evalHessian(f, x, varargin)
% Calculates the Hessian (i.e., d^2f/(dx_idx_j)) of the function f at point x
% according to Ridders' method:
%
% Ridders, CJF. (1982). Accurate computation of F'(x) and F'(x) F''(x). Advances in Engineering
%     Software, 4(2), 75-6.
%
% INPUT:
%    f             Function handle of a scalar real function of a real vector variable
%                  which are passed as *one* vector with *d* elements
%    x             Point at which to differentiate f
%
% OUTPUT:
%    Hessian       a set of 2nd derivatives of a given function f with
%                  respect to a real vector x.
%    min_err       Error estimate, a matrix with the same shape of Hessian
%
% OPTIONS:
%    Optionally, the third argument of the function can be a structure containing further
%    settings for Ridder's method.
%
%--------------------------------------------------------------------------
%    varargin{1}.init_eps       Initial  finite difference (default: 1)
%    varargin{1}.ratio          which is used to reduce eps on each step 
%                               (default: 0.8)
%    varargin{1}.min_steps      Minimum number of steps in eps (default: 5)
%    varargin{1}.max_steps      Maximum number of steps in eps (default: 100)
%    varargin{1}.threshold      Terminate if last step worse than preceding
%                               by a factor of threshold
%                               (default: 1.5)
%    varargin{1}.tolerance      terminal computation if the fact error is
%                               less than tolerance (default: 0.001)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    AUTHER: Changbo Zhu
%    E-mail: changbozhu@outlook.com           
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    d = length(x);
    I = eye(d);
    % results
    Hessian  = zeros(d,d);
    min_err  = realmax*ones(d,d);
    
    % Biuld struct of configuration for evaluating gradient of function f.    
    if nargin > 2
        options = varargin{1};
        options = evalHessian_options(options);
    else
        options = evalHessian_options();
    end
    if options.mpi
        % Loop through each component of variable x by using parfor     
        startmatlabpool;
        parfor i = 1:d
            Hessian_col = zeros(d,1);
            min_err_col = zeros(d,1);
            %ei = I(i,:);
            ei = zeros(1, d);
            ei(i) = 1;
            for j = 1:i
                %ej = I(:,j);
                ej = zeros(d, 1);
                ej(j) = 1;
                phi_eps = @(eps) phi(eps,f,x,ei,ej);
                [deriv2_ei_ej , min_err_ei_ej] = evalderiv2(phi_eps, options);
                Hessian_col(j) = deriv2_ei_ej;
                min_err_col(j) = min_err_ei_ej;
            end
            Hessian(:,i) = Hessian_col;
            min_err(:,i) = min_err_col;
        end
        Hessian = Hessian + Hessian' - diag(diag(Hessian));
        min_err = min_err + min_err' - diag(diag(min_err));
        closematlabpool;
    else
        for i = 1:d
            ei = I(i,:);
            for j = 1:i
                ej = I(:,j);
                phi_eps = @(eps) phi(eps,f,x,ei,ej);
                [deriv2_ei_ej , min_err_ei_ej] = evalderiv2(phi_eps, options);
                Hessian(i,j) = deriv2_ei_ej;
                min_err(i,j) = min_err_ei_ej;
                if j<i
                    Hessian(j,i) = deriv2_ei_ej;
                    min_err(j,i) = min_err_ei_ej;
                end
            end
        end         
    end
    
 end
 
 
 
 
 function [deriv2 , min_err] = evalderiv2(phi, options)
    
    deriv2 = NaN;
    
    % upacking parameters    
    min_err   = options.min_err;
    tolerance = options.tolerance;
    init_eps  = options.init_eps;
        
    ratio = options.ratio;       
    min_steps = options.min_steps;        
    max_steps = options.max_steps;
                
    threshold = options.threshold;
    
    % Initialize recording matrix A by max_steps times max_steps
    A = NaN(max_steps);
    % |------|--------------|---------------|---------------|--------------
    % |      |O(eps^(2*0+2))|O(eps^(2*1+2)) |O(eps^(2*2+2)) |O(eps^(2*n+2))
    % |------|--------------|---------------|---------------|--------------
    % | m=0: |  A0(eps)     |               |               |
    % | m=1: |  A0(r*eps)   |   A1(eps)     |               |
    % | m=2: |  A0(r^2*eps) |   A1(r*eps)   |   A2(eps)     |
    % | m=3: |  A0(r^3*eps) |   A1(r^2*eps) |   A2(r*eps)   |
    % |  :   |      :       |          :    |       :       |
    % |  :   |      :       |          :    |       :       |
    % | m=M; | A0(r^M*eps)  |A1(r^(M-1)*eps)|A2(r^(M-2)*eps)|
    % |------|--------------|---------------|---------------|--------------
    %  N = M = max_steps - 1  
    %  A(m,0) = (phi(r^m *eps) - phi(- r^m * eps))/(2*r^m *eps)
    %  A(m,n) = (A(m,n-1) -  coef*A(m-1 ,n-1))/(1 - coef), 0 < n < m
    %  absolutely_error(i,j) 
    %  = max(abs(A(m,n) -  A( m,n-1)), abs(A(m,n) -  A(m-1,n-1)))/(1 - r^2)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Initialize finite difference step epsilon (eps)
    eps = init_eps;
    M = max_steps - 1;
    N = max_steps - 1;
    % Compute central difference formula at initial step
    A(1,1) = phi(eps);
    
    % Loop to reduce eps 
    for m = 1:M
        n = 0;
        % Generate a new step size eps_m = ratio^m * init_eps
        eps = ratio*eps;
        
        % Compute 2nd order approximation by central difference formula at 
        % this step 
        A(m+1,n+1) = phi(eps);
        
        % Use square of ratio for extrapolation because errors increase
        % quadratically with eps (here, of course, they decrease quadratically
        % because we're reducing eps ...)
        ratio_sq = ratio^2;
        coeff = ratio_sq;
        % Fill the current row using Richardson extrapolation
        for n = 1:m            
            % Richardson extrapolation
            A(m+1,n+1) = ( A(m+1,n) - coeff*A(m,n) )/(1 - coeff);
            
            % Increment extrapolation factor
            coeff = coeff*ratio_sq;
            
            % Error on this trial is defined as the maximum absolute difference
            % to the extrapolation parents
            err_mn = max(abs(A(m+1,n+1)-A(m+1,n)),abs(A(m+1,n+1) - A(m,n)));
            
            if err_mn < min_err
                min_err = err_mn;
                deriv2 = A(m+1,n+1);
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
 
 function y = phi(eps,f,x,ei,ej)
    d = length(x);
    x = reshape(x,d,1);
    ei = reshape(ei,d,1);
    ej = reshape(ej,d,1);
    y = f(x + eps*ei + eps*ej ) - f(x - eps*ei + eps*ej ) ...
        - f(x + eps*ei - eps*ej ) + f(x - eps*ei - eps*ej );
    y = y/(4*eps^2);
 end
