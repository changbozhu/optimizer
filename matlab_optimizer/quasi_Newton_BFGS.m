function [ argMinVal, optims] = quasi_Newton_BFGS(obj, opt_vars0, varargin)
% This function implements the quasi-Newton minimization algorithm based on
% Line Search Method with Wolfe Condition.
% This variant of Newton Method is introduced by Broyden,Fletcher,Goldfarb,
% and Shanno, named as quasi Newton BFGS Algorithm.
%
% REFERENCE:
% [1]. Numerical Optimization (second edition), Jordge Nocedal, Stephen J.
%      Wright, Springer, New York, 2006
%--------------------------------------------------------------------------
% INPUT:
%    obj          A struct argument 
%        obj.obj_func     Function handle of the objective function to be optimized
%        obj.obj_grad     OPTIONAL FIELD.
%                         Function handle of the gradient of the objective function to be optimized
%    opt_vars0    The point at which the algorithm starts to search
%                 optimal solution.
%     varargin    Optional settings structure that can contain the
%                 following fields:
%       tolGrad   Convergence tolerance in the gradient
%       tolArg    Convergence tolerance in the argument
%       maxIter   Maximum number of iterations
%       maxRegu   Maximum number of regularizations
%       verbose   Boolean flag to turn output on (true) or off (false)
%
% OUTPUT:
%     argMinVal     Solution for the point of minimal objective function value
%     optims        Structure containing results in the following fields
%        minVal        The value of the function at its minimum
%        argMinVal     The argument of the function at its minimum
%        invHess       The inverse Hessian at the minimum calculated as a
%                      byproduct of optimization
%
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    AUTHER: Changbo Zhu
%    E-mail: changbozhu@outlook.com           
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nargin > 2
        options = varargin{1};
        evalGradOpts = options.evalGradOpts;
        options = quasi_Newton_BFGS_options(options);
    else
        options = quasi_Newton_BFGS_options();
        evalGradOpts = options.evalGradOpts;
    end   
    
    %quasi-Newton Method       
    % check gradient function and determin the parameter obj_para.objgrad_name
    if ~isfield(obj, 'obj_grad')
         obj.obj_grad = @(x) evalgradient(obj.obj_func, x, evalGradOpts);
    end
    
    if ~isfield(obj, 'obj_func')
         error('You must input an objective function via the struct argument obj.');
    end
    min_iters  = options.min_iters;
    max_iters  = options.max_iters;
    verbose    = options.verbose;
    step_mode  = options.step_mode;
    
    
    if verbose
        disp('  ')
        disp('Optimization Package is coded by Changbo Zhu.')
        disp(['-quasi_Newton_BFGS: '  datestr(now,'yyyy-mm-dd HH:MM:SS')])
    end
    
    % Make sure init is a column vector
    if ~iscolumn(opt_vars0)
        opt_vars0 = opt_vars0';
        if ~iscolumn(opt_vars0)
            error('Initial point has to be a row vector.');
        end
    end
    
    tstart = tic;
    if evalGradOpts.mpi
        startmatlabpool;
    end
    
    xk = opt_vars0;
    obj_fk  = obj.obj_func(xk);
    
    if verbose
        disp(' ')
        disp(['Initial argument: ', num2str(xk')])
        disp(['Initial value: ' num2str(obj_fk)])
    end
    
    grad = obj.obj_grad(xk);
    slope = - grad'*grad;
    
    dims = length(xk) ;
    I    = eye( dims );
    
 
    % BFGS Method Algorithm
        %     Given starting point x0, convergence tolerance epsilon >0,
        %           inverse Hessian approximation B0;
        %     iter = 0;
        %     While ||gradient(f(x(k)))|| > epsilon
        %           Compute search direction 
        %                                p(k) = - B(k) * gradient(f(x(k)));
        %           Set x(k+1) = x(k) + alpha(k) * p(k) where alpha_k is computed
        %           from a line search procedure to satisfy the (strong) 
        %           Wolfe conditions ;
        %
        %           Define s(k) = x(k+1)  -  x(k) ;
        %                 and y(k) = gradient(f(x(k+1))) - gradient(f(x(k)));
        %           Compute B(k+1) by means of BFGS update equation;
        %           iter = iter + 1;
        %      end(while)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-------------Start the quasi-Newton Method with BFGS--------------------
    tolgrad= options.tolgrad;
    tolarg= options.tolarg;  
        
    %---------Initialize inverse Hessian approximation B0-----------        
    % the initial matrix B0 often is set to some mutiple beta*I of the
    % identity,but there is no good general strategy for choosing the 
    % multiple beta.We ask the user to prescribe  a value delta for 
    % the norm of the first step,and then set 
    %                 B0 = delta / || g0 || * I 
    % to achieve this norm.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta = 1 ;               % Give the norm of the first step
    norm_grad = norm(grad);
    if norm_grad ==0
        norm_grad = 1;
    end
    beta  = delta/norm_grad;
    Bk = beta * I ;            % Generate the initial inverse Hessian 
                               % approximation
         
    % Do loop 
    %newton_iter = 0;
    converged_flag = 0; % 0--max_iter, 1-- tolgrad, 2-- tolarg
 
    for newton_iter = 1:max_iters  %(( norm(grad) > tolgrad )&&( newton_iter <= max_iters )) || ( newton_iter <= min_iters )
        if ( newton_iter > min_iters ) && max(...
            abs(grad')*max(abs(xk),1)./max(abs(obj_fk),1)...
            ) < tolgrad && max(abs(grad))<tolgrad
            if verbose
                disp(' ')
                disp('Converged on gradient size')
            end
            converged_flag = 1;
            break;
        end
        
        if newton_iter == max_iters
            disp('  ')
            disp('Please notice:')
            disp('Warning: optimization terminated because the maximum number of iterations is reached.')
            disp(['Warning: the norm of gradient at terminated point is ' num2str(norm(grad)) '>' num2str(tolgrad)  '!'])
            argMinVal = xk;
            minGrad   = grad;
        end
        % Compute the descent direction
        pk = - Bk * grad;
        %disp(['norm_pk='  num2str(norm(pk)) '; slope=' num2str(pk'*grad) ]
        if strcmpi(step_mode, 'regstep') 
            disp('This mode is developing...')
        elseif strcmpi(step_mode, 'linsearch') || strcmpi(step_mode, 'linesearch') 
            alpha =  Linesearch(obj.obj_func, xk, pk, obj.obj_grad, options );
        elseif strcmpi(step_mode, 'fix')
            alpha =1;
        end
        
        dx=alpha*pk;
        
        if verbose
            disp(['---quasi Newton BFGS Iteration:' num2str(newton_iter)])
            disp(['The norm of the current gradient ' num2str(norm(grad))])
            %disp(grad)
            disp(['Current step size: alpha=' num2str(alpha)])
            disp('dx=')
            disp(dx)
        end
        
        % Test for convergence
        if ( newton_iter > min_iters ) && (max(abs(dx)./abs(max(xk,1))) < tolarg)
            if verbose
                disp(' ')
                disp('Converged on step size')
            end
            converged_flag = 2;
            break;
        end
        
        xk = xk + dx;
        grad_old = grad;
        grad = obj.obj_grad(xk);
          
              
        % Define 
        %         sk = xk+1 - xk = alpha*pk 
        % and 
        %         yk = gradient(xk+1) - gradient(xk)
        % 
        sk = alpha*pk ;
        yk = grad - grad_old;
              
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % A heuristic that is often quite effective  is to scale the
        % starting matrix after the first step has been computed but
        % before the first BFGS update is performed. We change the 
        % provisional value B0 = I by setting
        %                    yk'sk
        %              B0 = ------- * I
        %                    yk'yk
        % before applying the update BFGS equation to obtain Bk.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if newton_iter == 0
             coff =( yk' * sk) / ( yk' * yk) ;
             Bk = coff * I ;
        end             
                          
        Bk = bfgs_eq(sk , yk , Bk);              
   end  % end(for-loop)
   
   
   
   if converged_flag ==1
        argMinVal = xk;
        minGrad   = grad;
        if verbose
            disp(['    -Converged at the ' num2str(newton_iter) 'th iteration.'])
            disp(['    -The norm of gradient at terminated point is ' num2str(max(abs(grad))) '<' num2str(tolgrad)  '.'])
        end
   elseif converged_flag ==2
       argMinVal = xk;
       minGrad   = grad;
       if verbose
            disp(['    -Converged at the ' num2str(newton_iter) 'th iteration.'])
            disp(['    -The difference of arguments at terminated point is ' num2str(max(abs(dx))) '<' num2str(tolarg)  '.'])
       end
   end
   telapse = toc(tstart);
   if verbose
        disp(['    -Optimization Over: ' datestr(now,'yyyy-mm-dd HH:MM:SS')])
        disp(['    -Elapsed time:' num2str(telapse) 's.'])
   end
   if evalGradOpts.mpi
       closematlabpool;
   end
   % Collect results
   optims = struct;
   optims.argMinVal = argMinVal;
   optims.minVal  =  obj.obj_func(argMinVal);
   optims.minGrad =  minGrad;
   optims.B_BFGS =  Bk;
   return;
end





function [ B ] = bfgs_eq( sk , yk , Bk )
% BFGS update equation:
%                        1
%                 rho = ------
%                      yk'sk
% and
%          Bk+1 = (I - rho*sk*yk')*Bk*(I - rho*yk*sk') + rho*sk*sk'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if yk'*sk > sqrt(eps*(yk'*yk)*(sk'*sk))
        dims = length( sk ) ;
        I    = eye( dims );
        
         rho = 1 / (yk'*sk);
        zvec = I - rho*sk*yk';
        B = zvec*Bk*zvec' + rho*(sk*sk');
     else
         B=Bk;
     end
     return;
end
