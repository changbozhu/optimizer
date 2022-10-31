function [ xBest, optims] = Particle_Swarm_Optimization(obj, x0, varargin)
% [xBest, fBest, info, dataLog] = Particle_Swarm_Optimization(objFun, x0, xdims, options)
%
% Particle Swarm Optimization
%
% This function minimizes or maximizes OBJFUN using some variants of particle 
% swarm optimization (e.g. PSO, EPSO, AEPSO). The optimization uses an initial 
% guess X0, and searches over a search space bounded by xLower and xUpper.
%
% INPUTS:
%   obj = objective function struct 
%         obj.obj_func  = objective function handle:
%                         F = obj.obj_func(x): the handle of the fitness function
%         obj.obj_grad  = gradient of objective function handle:
%       
%           x <--> [xdims, numP] : search point in n-dimensional space (for numP points)
%           F <--> [1, numP]     : objective function value, for each of m points
%   x0 <--> [xdims, 1] : initial search location
%                    Optional input. Set x0 = [] to ignore.
%   
%   options = option struct. All fields are optional, with defaults:
%       .xLower <--> [xdims, 1] : lower bounds on search space
%       .xUpper <--> [xdims, 1] : upper bounds on search space
%       .w_init        : search weight on current search direction (default: 0.4)
%       .c1_init       : search weight on global best (default: 0.9)
%       .c2_init       : search weight on local best (default: 0.9)
%       .numParticles  : population count  (default:numP= 3*xdims)
%       .maxIters      : maximum number of generations (default: 500)
%       .tolFun        : exit when variance in objective is < tolFun (default:1e-3)
%       .tolX          : exit when norm of variance in state < tolX (default:1e-3)
%       .flagVectorize : is the objective function vectorized? (default: false)
%       .flagMinimize = true = minimize objective
%           --> Set to false to maximize objective
%       .flagWarmStart = false = directly use initial guess?
%           --> true:  first particle starts at x0
%           --> false: all particles are randomly selected
%       .guessWeight : trade-off for initialization; range [0, 0.9) (default: 0.2)
%           --> 0.0  ignore x0; use random initialization [xLower, xUpper]
%           --> 0.9  heavy weight on initial guess (x0)
%       .display = 'iter';
%           --> 'iter' = print out info for each iteration
%           --> 'final' = print out some info on exit
%           --> 'off' = disable printing
%       .printMod = 1   (only used if display == 'iter')
%
% OUTPUTS:
%   xBest <--> [xdims, 1] : best point ever found
%   optims.argMinVal = xBest
%   optims.MinVal =  fBest <--> [1, 1]     : value of best point found
%   optims..exitFlag = how did optimization finish
%           0 = objective variance < tolFun
%           1 = reached max iteration count
%           2 = norm of state variance < tolX
%   optims.info = output struct with solver info
%       .input = copy of solver inputs:
%           .objFun
%           .x0
%           .xLower
%           .xUpper
%           .options
%       
%       .fEvalCount : how many calls to the objective function?
%       .X_Global     <--> [xdims,iter]  : best point in each generation
%       .F_Global     <--> [1,iter]      : value of the best point ever
%       .Pidx_Global  <--> [1,iter]      : index of the best point ever
%       .X_Best_Var   <--> [xdims,iter]  : variance in best point along each dim
%       .X_Var        <--> [xdims,iter]  : variance in current search along each dim
%       .X_Best_Mean  <--> [xdims,iter]  : mean in best point along each dim
%       .X_Mean       <--> [xdims,iter]  : mean in current search along each dim
%       .F_Best_Var   <--> [1,iter]      : variance in the best val at each gen
%       .F_Var        <--> [1,iter]      : variance in the current val at each gen
%       .F_Best_Mean  <--> [1,iter]      : mean of the population best value
%       .F_Mean       <--> [1,iter]      : mean of the current population value
%
%   dataLog(iter) = struct array with data from each iteration
%       .X          <--> [xdims,numP] : current position of each particle
%       .V          <--> [xdims,numP] : current "velocity" of each particle
%       .F          <--> [1,numP]     : value of each particle
%       .X_Best     <--> [xdims,numP] : best point for each particle
%       .F_Best     <--> [1,numP]     : value of the best point for each particle
%       .X_Global   <--> [xdims,1]    : best point ever (over all particles)
%       .F_Global   <--> [1,1]        : value of the best point ever
%       .I_Global   <--> [1,1]        : index of the best point ever
%
% NOTES:
%   This function uses a slightly different algorithm based on whether or
%   not the objective function is vectorized. If the objective is
%   vectorized, then the new global best point is only computed once per
%   iteration (generation). If the objective is not vectorized, then the
%   global best is updated after each particle is updated.
%
%
% REFERENCES:
%
%
% CHANGE LOG:
%
%  16th June, 2021
%   --> Add a function to evaluate a positive Hessian approximation by means
%       of BFGS method
%   ----------------------------------------------------------------------
%   15th Mary, 2021
%   --> Release the initial version
%

tic;
time_tag = strrep(datestr(now), ' ','');
time_tag = strrep(time_tag, ':','');

% check objective function and gradient function
if ~isfield(obj, 'obj_func')
     error('You must input an objective function via the struct argument obj.');
end


%%%% Options Struct:
xdims =length(x0);
if nargin > 2
    options = varargin{1};    
    options = Particle_Swarm_Optimization_Options(xdims, options);
else
    options = Particle_Swarm_Optimization_Options(xdim);
end

xLower = options.xLower;
xUpper = options.xUpper;
xdims  = options.xdims;

% Column vector validation:
[~, m] = size(x0);
if  m > 1
    error('x0 is not a valid size! Must be a column vector.')
end
[nRow, nCol] = size(xLower);
if nRow ~= xdims || nCol ~= 1
    error(['xLower is not a valid size! Must be [' num2str(xdims) ', 1]']);
end
[nRow, nCol] = size(xUpper);
if nRow ~= xdims || nCol ~= 1
    error(['xUpper is not a valid size! Must be [' num2str(xdims) ', 1]']);
end



%%%% Minimization problem or Maximization problem:
if options.flagMinimize
    optFun = @min;
else
    optFun = @max;
end

%%%% Check to see if user defined x0. If not, force defaults
if isempty(x0)
   x0 = 0.5*xLower + 0.5*xUpper;
   options.flagWarmStart = false;
end


%%%% Initialize the population

% Sample two random points in the search space for each particle
numP = options.numParticles;  % population size
XLower = xLower*ones(1,numP); % copy a column vector to a matrix
XUpper = xUpper*ones(1,numP); % copy a column vector to a matrix
% 2fish
%  X1 = -2 + 4.*rand(xdims,numP);
%  X2 = -2 + 4.*rand(xdims,numP);

%5fish
X1 = -1 + 2.*rand(xdims,numP);
X2 = -1 + 2.*rand(xdims,numP);

% Move initial points towards initial guess, by convex combination
% for initialization
X0 = x0*ones(1,numP); % copy a column vector to a matrix
X1 = options.guessWeight*X0 + (1 - options.guessWeight)*X1;
X2 = options.guessWeight*X0 + (1 - options.guessWeight)*X2;

% Initialize population:
X = X1;     % Initial position of the population
V = X2 - X1;  % Initial "velocity" of the population

% Check for warm start. If so, override random initial point with x0
if options.flagWarmStart
   X(:,1) = x0; 
   V(:,1) = zeros(size(x0));
end

if options.flagVectorize   % Batch process objective
    F = obj.obj_func(X);  % Function value at each particle in the population
else  % Objective not vectorized
    F = zeros(1,numP);
    for idx = 1:numP   % Loop over particles
        F(1,idx) = obj.obj_func(X(:,idx));
    end
end

X_Best = X;  % Best point, for each particle in the population
F_Best = F;  % Value of best point, for each particle in the population


[F_Global, Pidx_Global] = optFun(F_Best); % Value of best point ever, over all points
X_Global = X(:, Pidx_Global); % Best point ever, over all  points
% X_Global_Record = X_Global;
% update_iter = 1;

%%%% Allocate memory for the dataLog
maxIters = options.maxIters;
iter =options.maxRecIters;
dataLog(options.maxRecIters) = makeStruct(iter, X, V, F, X_Best, F_Best, X_Global, F_Global, Pidx_Global);
%%%% Allocate memory for info
info = makeInfo(options);


%%%% MAIN LOOP:
exitFlag = 1;   %Assume that we will reach maximum iteration
numdataLog = 1;
dataLog_idx = 1;
t_iter=0;

for iter = 1:maxIters
    t1=toc;
    %%% Compute new generation of points:
    if iter > 1   % Then do an update on each particle
        r1 = rand(xdims,numP);
        r2 = rand(xdims,numP);
        
        switch options.variant
            case 'PSO'
                w = weight_decay(iter,options);
                c1 = options.c1_init;
                c2 = options.c2_init;
                X_a = X_Best;
            case 'EPSO'
                w = weight_decay(iter,options);
                c1 = options.c1_init;
                c2 = options.c2_init;
                X_a = mean(X_Best,2)*ones(1, numP);
            case 'AEPSO'
                w = weight_decay(iter, options);
                [c1, c2] = c_reg(F, F_Global, options);
                X_a = mean(X_Best, 2)*ones(1, numP);
        end
        
        t2=toc;
        if options.flagVectorize   % Batch process objective            
            V =  ...   %Update equations
                w*V + ...    % Current search direction
                c2*r2.*((X_Global*ones(1,numP)) - X) + ...  % Global direction
                c1*r1.*(X_a - X);    % Local best direction
            X_New = X + V;  % Update position
            X = max(min(X_New, XUpper), XLower);   % Clamp position to bounds
            
            F = obj.obj_func(X);   %Evaluate
            
            F_Best_New = optFun(F_Best, F);   %Compute the best point
            idxUpdate = F_Best_New ~= F_Best;  % Which indicies updated?
            X_Best(:,idxUpdate) = X(:,idxUpdate);  %Copy over new best points
            F_Best = F_Best_New;
            [F_Global, Pidx_Global] = optFun(F_Best); % Value of best point ever, over all points
            X_Global = X(:, Pidx_Global); % Best point ever, over all  points
%             X_Global_Record = [X_Global_Record, X_Global];  
%             update_iter = iter;
            
        else   %Objective is not vectorized.            
            for idx = 1:numP   %%%%%%% Loop over particles  %%%%%%%%%%%%%%
                
                V(:,idx) =  ...   %Update equations
                    w*V(:,idx) + ...    % Current search direction
                    c2*r2(:,idx).*(X_Global-X(:,idx)) + ...  % Global direction
                    c1*r1(:,idx).*(X_a(:,idx)-X(:,idx));    % Local best direction
                X_New = X(:,idx) + V(:,idx);  % Update position
                X(:,idx) = max(min(X_New, xUpper), xLower);   % Clamp position to bounds
                
                F(:,idx) = obj.obj_func(X(:,idx));   %Evaluate
                
                [F_Best(1,idx), iBest] = optFun([F(1,idx), F_Best(1,idx)]);
                if iBest == 1  %Then new point is better!
                    X_Best(:,idx) = X(:,idx);
                    [F_Global, iBest] = optFun([F_Best(1,idx), F_Global]);
                    if iBest == 1 %Then new point is the global best!
                        X_Global = X_Best(:,idx);
%                         X_Global_Record = [X_Global_Record, X_Global]; 
%                         update_iter = iter;
                    end
                
                end
                
            end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
%         if iter - update_iter > 10 && norm(X_Global_Record(:,end) - X_Global)>sqrt(options.tolX)
%             X_Global_Record = [X_Global_Record, X_Global]; 
%             update_iter = iter;
%         end
        t3=toc;
        t_iter = 0.9*t_iter+ 0.1*(t3-t2) ;
        %disp(['    Evaluate elpsed time: ' num2str( (maxIters-iter)*t_iter/3600) 'hours.' 't3-t2=' num2str(t3-t2)])
    end
    
    %%% Log Data
    dataLog(dataLog_idx) = makeStruct(iter, X, V, F, X_Best, F_Best, X_Global, F_Global, Pidx_Global);
    info = logInfo(info, iter, X, F, X_Best, F_Best, X_Global, F_Global, Pidx_Global);
    
%     %%% Plot
%     if ~isempty(options.plotFun)
%         options.plotFun(dataLog(iter), iter);
%     end
    
    %%% Print:
    xVar = norm(info.X_Var(:,end));
    if strcmp('iter',options.display)
        if mod(iter-1,options.printMod)==0
            fprintf('iter: %3d, fBest: %9.3e, F_Mean: %9.3e,  fVar: %9.3e   xVar: %9.3e  remaining time: %3.2f hours.  \n',...
                iter,  info.F_Global(end), info.F_Mean(end), info.F_Var(1,end), xVar, (maxIters-iter)*t_iter/3600);
        end
    end
    
    %%% Convergence:
    if info.F_Var(1,end) < options.tolFun
        exitFlag = 0;
        %info = truncateInfo(info,maxIters,iter);
        break
    elseif xVar < options.tolX
        exitFlag = 2;
        %info = truncateInfo(info,maxIters,iter);
        break
    end
    
    if dataLog_idx >= options.maxRecIters || iter >= maxIters
        dataLog = dataLog(1:dataLog_idx);
        mkdir(['./dataLog/' time_tag ]);
        save(['./dataLog/' time_tag '/dataLog' num2str(numdataLog) '.mat'],'dataLog')
        numdataLog = numdataLog + 1;
        dataLog_idx = 1;
    else
        dataLog_idx = dataLog_idx + 1;
    end
end


xBest = info.X_Global(:,end);
fBest = info.F_Global(end);
info.input = makeStruct(obj, x0, xdims, options);  %Copy inputs
info.fEvalCount = iter*numP;

% compute Hessian approximation by means of BFGS Equation
% if options.HessianApproxMod
%     if  ~isfield(obj, 'obj_grad')
%         evalGradOpts = struct('mpi', true, ...
%                           'init_eps', 0.01, ...
%                           'ratio', 0.8, ...
%                           'min_steps', 5, ...
%                           'max_steps', 100, ...
%                           'threshold', 1.5, ...
%                           'tolerance', 1e-3 ...
%         );
%         obj.obj_grad = @(x) evalgradient(obj.obj_func, x, evalGradOpts);
%         if evalGradOpts.mpi
%             startmatlabpool;
%         end
%     end
% end
% Grad_Global = obj.obj_grad(X_Global_Record(:,1));
% B_BFGS=10/norm(Grad_Global)*eye(xdims);
% K =  size(X_Global_Record, 2);
% for i = 2:1:K
%     Grad_Global_old = Grad_Global;
%     Grad_Global = obj.obj_grad(X_Global_Record(:,i));
%     sk = X_Global_Record(:,i) - X_Global_Record(:,i-1);
%     yk =   Grad_Global - Grad_Global_old;
%     B_BFGS = B_BFGS  ...
%          - (B_BFGS * sk * sk' * B_BFGS)/(sk' * B_BFGS * sk) ...
%          + yk*yk' / (yk' * sk);
% end
% if isfield(evalGradOpts,'mpi') && evalGradOpts.mpi
%     closematlabpool;
% end


% create a struct optims
optims = struct();
optims.argMinVal =  xBest;
optims.minVal    =  fBest;
optims.info      =  info;
optims.dataLog   =  dataLog;
optims.exitFlag  =  exitFlag;
%optims.B_BFGS = B_BFGS;
%%% Print:
if strcmp('iter',options.display) || strcmp('final',options.display)
    switch exitFlag
        case 0
            fprintf('Optimization Converged. Exit: fVar < tolFun \n');
        case 1
            fprintf('Maximum iteration reached. \n');
        case 2
            fprintf('Optimization Converged. Exit: norm(xVar) < tolX \n');
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = weight_decay(iter,options)
w_max  = options.w_init;
w_min  = options.w_min;
maxIters = options.maxIters;
switch options.w_decay_model
    case 'linear'
        w = w_max-  iter/maxIters*( w_max - w_min );
    case 'quadratic'
        w = w_max-  (iter/maxIters)^2*( w_max - w_min );
end
return; 
end
function [c1, c2] = c_reg(F, F_Global,options)
c_alpha = options.c_alpha;
Fa = mean(F);
ck =exp( - c_alpha*abs(Fa - F_Global) );
c1 = 1 - ck;
c2 = ck;
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function info = makeInfo(options)
%
% Initializes a recorder struct info
xdims = options.xdims;
maxIters = options.maxIters;
maxRecIters = options.maxRecIters;
maxSteps = min([maxRecIters, maxIters]);
info.X_Global = zeros(xdims,maxSteps);
info.F_Global = zeros(1,maxSteps);
info.Pidx_Global = zeros(1,maxSteps);
info.X_Best_Var = zeros(xdims,maxSteps);
info.F_Best_Var = zeros(1,maxSteps);
info.X_Best_Mean = zeros(xdims,maxSteps);
info.F_Best_Mean = zeros(1,maxSteps);
info.X_Var = zeros(xdims,maxSteps);
info.F_Var = zeros(1,maxSteps);
info.X_Mean = zeros(xdims,maxSteps);
info.F_Mean = zeros(1,maxSteps);
info.iter = NaN(1,maxSteps);
end
function info = logInfo(info, iter, X, F, X_Best, F_Best, X_Global, F_Global, Pidx_Global)
%
% Initializes a recorder struct info
info.X_Global(:,1) = [];
info.X_Global=[info.X_Global,  X_Global];
    
info.F_Global(:,1) = [];
info.F_Global=[info.F_Global,  F_Global];
        
info.Pidx_Global(:,1) = [];
info.Pidx_Global=[info.Pidx_Global,  Pidx_Global];
    
info.X_Var(:,1) = [];
info.X_Var = [info.X_Var, var(X,0,2)];
    
info.X_Best_Var(:,1) = [];
info.X_Best_Var = [info.X_Best_Var, var(X_Best,0,2)];
    
info.X_Mean(:,1) = [];
info.X_Mean = [info.X_Mean, mean(X, 2)];
    
info.X_Best_Mean(:,1) = [];
info.X_Best_Mean = [info.X_Best_Mean, mean(X_Best, 2)];
    
info.F_Var(:,1) = [];
info.F_Var = [info.F_Var, var(F)];
    
info.F_Best_Var(:,1) = [];
info.F_Best_Var = [info.F_Best_Var, var(F_Best)];
    
info.F_Mean(:,1) = [];
info.F_Mean = [info.F_Mean, mean(F)];
    
info.F_Best_Mean(:,1) = [];
info.F_Best_Mean = [info.F_Best_Mean, mean(F_Best)];
    
info.iter(:,1) = [];
info.iter = [info.iter, iter];
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function info = truncateInfo(info,maxIters,iter)
% %
% % Removes the empty entries in the info struct
% 
% names = fieldnames(info);
% for i=1:length(names)
%     if (isnumeric(info.(names{i})))   % Check if it's a matrix
%         if size(info.(names{i}),2) == maxIters    % Check if it is iteration data
%             info.(names{i}) = info.(names{i})(:,1:iter);
%         end
%     end
% end
% 
% end




