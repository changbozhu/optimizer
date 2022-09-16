% TEST  --  PSO  -- Particle Swarm Optimization
%
% Test 4:  Styblinski-Tang function (2D)
%
% 

clc; clear; clear global; 
%figure(400); clf;

%%%% Set up problem

objFun = @StyblinskiTang;   % Minimize this function
x0 = [0;0];  % initial guess

options.w_init  = 0.5;  % weight on current search direction
options.c1_init = 1.1;   % weight on local best search direction
options.c2_init = 1.1;  % weight on global best search direction

options.variant= 'PSO';
options.numParticles = 50;
options.maxIters = 1000;
options.xLower = -5*ones(2,1); % lower bound on the search space
options.xUpper =  5*ones(2,1); % upper bound on the search space

%options.plotFun = @plotStyblinskiTang;  % Plots progress


%%%% Solve
[xBest, optims] = Particle_Swarm_Optimization(objFun, x0, options);

%%%% Analysis
figure(401); clf;
plotPsoHistory(optims.info);


