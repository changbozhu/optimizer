function [ options_default ] = Linesearch_options( varargin )
    % Configure parameters for Linesearch with Wolfe Conidtion
    % Please reference to the file Linesearch.m for More details 
    options_default = struct;
   
    options_default.c1 = 0.0001;
    options_default.c2 = 0.9;
    % Initialize minimal step length alpha0
    % or define the left point of the step length interval
    options_default.alpha0   = 0;
    % Initialize step length alpha_i which chose from (alpha0, alpha_max)
    options_default.alpha_i  = 1;
    % Initialize maximal step length alpha_max
    % or define the right point of the step length interval
    options_default.alpha_max = 120;
    % interpolation method in a search interval
    options_default.interp_method = 'bisection'; 
        
    if nargin > 0
        options = varargin{1};

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
    end
    return;
end

