% TEST  --  PSO  -- Particle Swarm Optimization
%
% Test 5:  Styblinski-Tang function (5-D)
%
% There are many local minimums to this problem, but only one global
% minimum. All are of similar value. 
%   >> help StyblinskiTang   % For more details
%

clc; clear; 

%%%% Set up problem

objFun = @StyblinskiTang;   % Minimize this function
x0 = -2*ones(5,1);  % initial guess

options.w_init  = 0.5;  % weight on current search direction
options.c1_init = 1.0;   % weight on local best search direction
options.c2_init = 1.0;  % weight on global best search direction

options.variant= 'PSO';
options.numParticles = 50;
options.maxIters = 1000;
options.xLower = -5*ones(5,1); % lower bound on the search space
options.xUpper =  5*ones(5,1); % upper bound on the search space

options.flagVectorize = true;  % Use vectorized objective function
options.flagWarmStart = true;  % Include x0 in first generation

%%%% Solve
[xBest, optims] = Particle_Swarm_Optimization(objFun, x0, options);

%%%% Analysis
figure(501); clf;
plotPsoHistory(optims.info);


