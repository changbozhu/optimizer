function [ options_default ] = quasi_Newton_BFGS_options_rosenbrock(varargin)
    %initialize data sructure for Line Search
    options_default = struct;
    
    
    % norm precision or minimal norm of the gradient, 
    % it should be close to 0.
    options_default.tolgrad = 1e-4;   
    % max iteration steps
    options_default.min_iters = 10;
    options_default.max_iters = 2000;
    options_default.verbose = true;
    % you must choose a string from 'linsearch' , 'regStep' or 'fix'
    options_default.step_mode = 'regStep'; 
   
    % Configure parameters for Linesearch and Wolfe Conidtion
    % Please reference to the file Linesearch.m for More details 
    options_default.c1 = 0.0001; % this argument should be close to 0; and its' upper bound is c2.
    options_default.c2 = 0.9;    % c2 is chosen from (c1, 1) and should be close to 1.
    
    % Initialize minimal step length alpha0
    % or define the left point of the step length interval
    options_default.alpha0   = 0; % It's a fixed parameter, please don't adjust it
    
    % Initialize step length alpha_i which chose from (alpha0, alpha_max)
    options_default.alpha_i  = 1;
    
    % Initialize maximal step length alpha_max
    % In Linesearch model, it defines the right point of the step length interval
    % In regularizing step length model, it defines the maximum step length
    % moving in the descent direction
    options_default.alpha_max = 120; % The default value is 120, you can adjust it.
    
    % Interpolation method in a search interval
    % The default method is 'bisection interpolation', you can add your 
    % optional method. Please reference to [1].
    options_default.interp_method = 'bisection';    
   
    %
    if nargin > 0
        options = varargin{1};
        if isfield(options,'min_iters')
            options_default.min_iters = options.min_iters;
        end
        if isfield(options,'max_iters')
            options_default.max_iters = options.max_iters;
        end
        
        if isfield(options,'step_mode')
            options_default.step_mode = options.step_mode;
        end
        
        if isfield(options,'tolgrad')
            options_default.tolgrad = options.tolgrad;
        end

        if isfield(options,'c1')
            options_default.c1 = options.c1;
        end
        
        if isfield(options,'c2')
            options_default.c2 = options.c2;
        end
        
        if isfield(options,'alpha0')
            options_default.alpha0 = options.alpha0;
        end
        
        if isfield(options,'alpha_i')
            options_default.alpha_i = options.alpha_i;
        end
        
        if isfield(options,'alpha_max')
            options_default.alpha_max = options.alpha_max;
        end
        
        if isfield(options,'interp_method')
            options_default.interp_method = options.interp_method;
        end
        
        if isfield(options,'verbose')
            options_default.verbose = options.verbose;
        end
    end
    
    evalGradOpts_default =struct;
    % Defaults
    evalGradOpts_default.init_eps        = 1;
    evalGradOpts_default.ratio           = 0.8;
    evalGradOpts_default.min_steps       = 5;
    evalGradOpts_default.max_steps       = 100;

    evalGradOpts_default.threshold       = 1.5;    
    evalGradOpts_default.min_err         = realmax;
    evalGradOpts_default.tolerance       = 1e-5;
    evalGradOpts_default.mpi             = false;
    
    % Overrides
    if nargin > 0
        options = varargin{1};
        if isfield(options,'evalGradOpts')
            evalGradOpts = options.evalGradOpts;
            if isfield(evalGradOpts,'mpi')
                options_default.mpi = evalGradOpts.mpi;
            end
            
            if isfield(evalGradOpts,'init_eps')
                evalGradOpts_default.init_eps = evalGradOpts.init_eps;
            end
            
            if isfield(evalGradOpts,'ratio')
                evalGradOpts_default.ratio = evalGradOpts.ratio;
            end
            
            if isfield(evalGradOpts,'min_steps')
                evalGradOpts_default.min_steps = evalGradOpts.min_steps;
            end
            
            if isfield(evalGradOpts,'max_steps')
                evalGradOpts_default.max_steps = evalGradOpts.max_steps;
            end
            
            if isfield(evalGradOpts,'threshold')
                evalGradOpts_default.threshold = evalGradOpts.threshold;
            end
            
            if isfield(evalGradOpts,'tolerance')
                evalGradOpts_default.tolerance = evalGradOpts.tolerance;
            end
        end
    end
    if evalGradOpts_default.ratio >= 1 || evalGradOpts_default.ratio <0
        error('Reduction ratio must be in (0,1).')
    end
    options_default.evalGradOpts = evalGradOpts_default;
end

