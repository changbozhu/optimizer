% TEST  --  PSO  -- Particle Swarm Optimization
%
% Test 2:  Himmelblau's function 
%
% 

clc; clear; clear global; 
%figure(200); clf;
global HimmelblauSaveAnimation    %Set to true to save the animation

%%%% Set up problem

objFun = @Himmelblau;   % Minimize this function
x0 = [0;0];  % initial guess

% options.w_init  = 0.4;  % weight on current search direction
% options.c1_init = 0.9;   % weight on local best search direction
% options.c2_init = 0.9;  % weight on global best search direction

options.variant= 'PSO';
options.numParticles = 50;
options.maxIters = 1000;
options.xLower = -1e75*ones(2,1); % lower bound on the search space
options.xUpper =  1e75*ones(2,1); % upper bound on the search space

%options.plotFun = @plotHimmelblau;  % Plots progress
HimmelblauSaveAnimation = false;   %Save an animation


%%%% Solve
[xBest, optims] = Particle_Swarm_Optimization(objFun, x0, options);

%%%% Analysis
figure(201); clf;
plotPsoHistory(optims.info);


