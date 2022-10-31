% TEST  --  PSO  -- Particle Swarm Optimization
%
% First test, simple quadratic bowl

clc; clear;

%%%% Set up problem

objFun = @(x)( sum(x.^2,1) );   % Minimize this function

x0 = [0;0];  % No initial guess
%options.xdims = 2;
options.w_init  = 0.4;  % weight on current search direction
options.w_min  = 0.01;  % weight on current search direction
% options.c1_init = 0.9;   % weight on local best search direction
% options.c2_init = 0.9;  % weight on global best search direction
options.c_alpha = 0.5;


options.variant= 'AEPSO';
options.numParticles = 50;
options.maxIters = 1000;
%options.xLower = -sqrt(realmax)/5* ones(options.xdims,1); % lower bound on the search space
%options.xUpper =  sqrt(realmax)/5* ones(options.xdims,1); % upper bound on the search space

%options.plotFun = @plotBowl;  % Plots progress


%%%% Solve
[xBest, optims] = Particle_Swarm_Optimization(objFun, x0, options);

%%%% Analysis
figure(101); clf;
plotPsoHistory(optims.info);


