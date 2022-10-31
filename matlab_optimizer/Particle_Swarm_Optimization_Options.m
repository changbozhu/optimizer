function options_default = Particle_Swarm_Optimization_Options(xdims, varargin)
%
%
% OPTIONS:
%    Optionally, the third argument of the function can be a structure containing further
%    settings for Ridder's method.
%    varargin{1}.xdims
%    varargin{1}.xLower          the lower bound on the position x (default: -reamax)
%    varargin{1}.xUpper          the upper bound on the position x (default: reamax)
%    varargin{1}.maxIters        maximum number of iterations in PSO (default: 500)
%    varargin{1}.tolFun          exit when variance in objective is < tolFun  (default: 0.001)
%    varargin{1}.tolX            exit when variance in position is < tolX  (default: 0.001)
%    varargin{1}.w_init 
%    varargin{1}.w_min
%    varargin{1}.c1_init        
%    varargin{1}.c2_init
%    varargin{1}.c_alpha
%    varargin{1}.numParticles   
%    varargin{1}.variant        
%    varargin{1}.flagWarmStart  directly use initial guess?
% --------------------------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    AUTHER: Changbo Zhu
%    E-mail: changbozhu@outlook.com           
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xLower = -1e6;
    xUpper = 1e6;
    options_default = struct(...
                        'xdims', xdims, ... 
                        'xLower',    xLower * ones(xdims,1),  ... % Lower bound of position. It must be a column vector.
                        'xUpper',    xUpper * ones(xdims,1),  ... % Upper bound of postion. It must be a column vector.
                        'maxIters', 500,  ... % maximum iterations
                        'maxRecIters', 50, ...
                        'tolFun', 1e-3,   ... % exit when variance in objective is < tolFun
                        'tolX', 1e-9,     ... % exit when variance in position is < tolX
                        'w_init',  0.9,   ... % maximumum and initial weight on current search direction
                        'w_min',  0.4,    ... % minimum weight on current search direction
                        'w_decay_model', 'quadratic', ...
                        'c1_init', 0.9,   ... % weight on local best search direction
                        'c2_init', 0.9,   ... % weight on global best search direction
                        'c_alpha', 0.5,   ...
                        'numParticles', 3*xdims,... % the number of particles
                        'variant', 'PSO',       ... % 'PSO', 'EPSO','AEPSO'
                        'flagMinimize', true,   ... % minimize objective, if set to false to maximize objective
                        'flagWarmStart',true,   ... % directly use initial guess?
                                                ... % ---true:  first particle starts at x0
                        'flagVectorize',false,  ... % Is the objective function vectorized?
                                                ... % ---false: all particles are randomly selected
                        'guessWeight', 0.2,      ...% trade-off for initialization; range [0, 0.9)
                                                ... % ---0.0  ignore x0; use random initialization [xLow, xUpp]
                                                ... % ---0.9  heavy weight on initial guess (x0)     
                        'display', 'iter',      ... % --- 'iter': print out info for each iteration  
                                                ... % --- 'final': print out some info on exit
                                                ... % ---'off' = disable printing
                        'printMod', 1,          ... % only used if display == 'iter';
                        'HessianApproxMod', 1   ...
                        );
    
    % Overrides    
    if nargin > 1
        options = varargin{1};        
        if xdims <= 0 && isfield(options,'xdims')
            options_default.xdims  =  options.xdims;
            options_default.xLower =   xLower * ones(options_default.xdims,1);
            options_default.xUpper =   xUpper * ones(options_default.xdims,1);
            options_default.numParticles = 3*options.xdims;
        end

        if isfield(options,'xLower')
            options_default.xLower = options.xLower;
        end
        
        if isfield(options,'xUpper')
            options_default.xUpper = options.xUpper;
        end
        
        if isfield(options,'maxIters')
            options_default.maxIters = options.maxIters;
        end
        
        if isfield(options,'maxRecIters')
            options_default.maxRecIters = options.maxRecIters;
        end
        
        if isfield(options,'tolFun')
            options_default.tolFun = options.tolFun;
        end
        
        if isfield(options,'tolX')
            options_default.tolX = options.tolX;
        end
        
        if isfield(options,'w_init')
            options_default.w_init = options.w_init;
        end
        
        if isfield(options,'w_min')
            options_default.w_min = options.w_min;
        end
        if isfield(options,'w_decay_model')
            options_default.w_decay_model = options.w_decay_model;
        end
        if isfield(options,'c1_init')
            options_default.c1_init = options.c1_init;
        end
        
        if isfield(options,'c2_init')
            options_default.c2_init = options.c2_init;
        end
        
        if isfield(options,'c_alpha')
            options_default.c_alpha = options.c_alpha;
        end
        
        if isfield(options,'numParticles')
            options_default.numParticles = options.numParticles;
        end
        
        if isfield(options,'variant')
            options_default.variant = options.variant;
        end
        
        if isfield(options,'flagMinimize')
            options_default.flagMinimize = options.flagMinimize;
        end
        
        if isfield(options,'flagWarmStart')
            options_default.flagWarmStart = options.flagWarmStart;
        end
        
        if isfield(options,'guessWeight')
            options_default.guessWeight = options.guessWeight;
        end
        
        if isfield(options,'display')
            options_default.display = options.display;
        end
        
        if isfield(options,'printMod')
            options_default.printMod = options.printMod;
        end
        
        if isfield(options,'HessianApproxMod')
            options_default.HessianApproxMod = options.HessianApproxMod;
        end
    end
    
    if options_default.xdims <= 0
        error('You must give a proper dimesion value of optimised variable or an initial variable.')
    end
end


