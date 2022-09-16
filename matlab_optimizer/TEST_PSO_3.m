% TEST  --  PSO  -- Particle Swarm Optimization
%
% Test 3:  Goldstein-Price function 
%
% 

clc; clear; clear global; 
%figure(300); clf;

%%%% Set up problem

objFun = @GoldsteinPrice;   % Minimize this function

x0 = [ ];  % initial guess

% options.w_init  = 0.4;  % weight on current search direction
% options.c1_init = 0.9;   % weight on local best search direction
% options.c2_init = 0.9;  % weight on global best search direction

% options.tolX = 1e-8;
% options.tolFun = 1e-4;

% options.flagVectorize = true;  % Objective function is vectorized

options.variant= 'PSO';
options.xdims = 2;
options.numParticles = 50;
options.maxIters = 1000;
options.xLower = -2*ones(2,1); % lower bound on the search space
options.xUpper =  2*ones(2,1); % upper bound on the search space

%options.plotFun = @plotGoldsteinPrice;  % Plots progress

%%%% Solve
[xBest, optims] = Particle_Swarm_Optimization(objFun, x0, options);

%%%% Analysis
figure(301); clf;
plotPsoHistory(optims.info);


