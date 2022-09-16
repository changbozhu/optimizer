function options_default = evalHessian_options(varargin)
% Evaluates derivative of the function f at point x according to Ridders' method:
%
% Ridders, CJF. (1982). Accurate computation of F'(x) and F'(x) F''(x). Advances in Engineering
%     Software, 4(2), 75-6.
%
% INPUT:
%    phi           Function handle of a scalar real function of one real variable
%
% OUTPUT:
%    deriv         First derivative of f at x
%    min_error     Error estimate
%
% OPTIONS:
%    Optionally, the third argument of the function can be a structure containing further
%    settings for Ridder's method.
%
%    varargin{1}.init_eps    Initial  finite difference (default: 1)
%    varargin{1}.ratio       which is used to reduce h on each step (default: 0.8)
%    varargin{1}.min_steps   Minimum number of steps in eps (default: 5)
%    varargin{1}.max_steps   Maximum number of steps in eps (default: 100)
%    varargin{1}.threshold   Terminate if last step worse than preceding 
%                            by a factor of threshold (default: 1.5)
%
% --------------------------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    AUTHER: Changbo Zhu
%    E-mail: changbozhu@outlook.com           
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options_default = struct;
    % Defaults
    options_default.init_eps        = 1;
    options_default.ratio           = 0.8;
    options_default.min_steps       = 5;
    options_default.max_steps       = 100;

    options_default.threshold       = 1.5;    
    options_default.min_err         = realmax;
    options_default.tolerance       = 1e-5;
    options_default.mpi             = false;
    
    % Overrides
    if nargin > 0
        options = varargin{1};
        if isfield(options,'mpi')
            options_default.mpi = options.mpi;
        end

        if isfield(options,'init_eps')
            options_default.init_eps = options.init_eps;
        end
        
        if isfield(options,'ratio')
            options_default.ratio = options.ratio;
        end
        
        if isfield(options,'min_steps')
            options_default.min_steps = options.min_steps;
        end
        
        if isfield(options,'max_steps')
            options_default.max_steps = options.max_steps;
        end
        
        if isfield(options,'threshold')
            options_default.threshold = options.threshold;
        end
        
        if isfield(options,'tolerance')
            options_default.tolerance = options.tolerance;
        end
    end
    if options_default.ratio >= 1 || options_default.ratio <0
        error('Reduction ratio must be in (0,1).')
    end
end


