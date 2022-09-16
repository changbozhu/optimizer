% TEST  --  PSO  -- Particle Swarm Optimization
%
% Test 6:  Styblinski-Tang function (5-D)
%
% There are many local minimums to this problem, but only one global
% minimum. All are of similar value.
%   >> help StyblinskiTang   % For more details
%
% This script is a meta-optimization, running a simple grid search to
% determine which set of paramters will cause the optimization to most
% rapidly converge to a local minimum.
%
% Note - the resulting set of parameters are not good in general, since
% they will always find the closest local minimum.
%
% This script will take some time to execute, since it runs many
% optimizations one after the other.
%

clc; clear;

%%%% Set up problem

objFun = @StyblinskiTang;   % Minimize this function
x0 = zeros(5,1);  % initial guess

options.variant= 'PSO';
options.maxIters = 100;
options.tolFun = 1e-6;
options.tolX = 1e-6;
options.flagVectorize = false;
options.guessWeight = 0.11;
options.display = 'off';
options.xLower = -5*ones(5,1); % lower bound on the search space
options.xUpper =  5*ones(5,1); % upper bound on the search space

%%%% Select the grid search parameters:
Z_w = linspace(0.1, 0.7, 9);
Z_c1c2 = linspace(0.7, 1.6, 9);
Z_numP = [8, 12, 16, 20, 24];

N_REPEAT = 5;  %Run optimization this many times for each set of params

nw = length(Z_w);
nc1c2 = length(Z_c1c2);
nPartIter = length(Z_numP);

nTrial = nw*nc1c2*nPartIter;
W = zeros(nTrial,1);
C1C2 = zeros(nTrial,1);
Part = zeros(nTrial,1);

F_Eval_Count = zeros(nTrial,1);
F_Best = zeros(nTrial,1);


nEval = zeros(1,N_REPEAT);
fVal = zeros(1,N_REPEAT);

idx = 0;
for i=1:nw
    for j=1:nc1c2
        for k=1:nPartIter
            idx = idx + 1;
            
            %%%% Unpack parameters:
            options.w_init = Z_w(i);
            options.c1_init = Z_c1c2(j);
            options.c2_init = Z_c1c2(j);
            options.numParticles = Z_numP(k);
            
            %%%% Log parameters (lazy way)
            W(idx) = Z_w(i);
            C1C2(idx) = Z_c1c2(j);
            Part(idx) = Z_numP(k);
            
            %%%% Solve
            for rr = 1:N_REPEAT
                [~, optims] = Particle_Swarm_Optimization(objFun, x0, options);
                nEval(rr) = optims.info.fEvalCount;
                fVal(rr) = optims.minVal;
            end
            F_Eval_Count(idx) = mean(nEval(rr));
            F_Best(idx) = mean(fVal);
            
            %%%% User Read-Out:
            fprintf('Iter:  %d / %d \n', idx, nTrial);
            
        end
    end
end

%%%% Data Analysis

FF = [F_Eval_Count, F_Best];

[~,IDX] = sortrows(FF,[-1,-2]);   %[worst --> best]

% for i=1:nTrial
%     fprintf('nEval: %4d,  fVal: %6.3e,  W: %4.2f, C1:  %4.2f,  nPart: %d  \n',...
%         F_Eval_Count(IDX(i)),  F_Best(IDX(i)),  W(IDX(i)), C1C2(IDX(i)), Part(IDX(i)));
% end


%%%% Agregate the top N parameter runs:
N = 10; ii = length(IDX) + 1 - (1:N);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i=1:N
    fprintf('nEval: %4d,  fVal: %6.3e,  W: %4.2f, C1:  %4.2f,  nPart: %d  \n',...
        F_Eval_Count(IDX(ii(i))),  F_Best(IDX(ii(i))),  W(IDX(ii(i))), C1C2(IDX(ii(i))), Part(IDX(ii(i))));
end
fprintf('nEval:  mean = %6.1f,  median = %6.1f, range = %6.1f \n', mean(F_Eval_Count(IDX(ii))), median(F_Eval_Count(IDX(ii))),range(F_Eval_Count(IDX(ii))));
fprintf('fVal:   mean = %6.1f,  median = %6.1f, range = %6.1f \n', mean(F_Best(IDX(ii))), median(F_Best(IDX(ii))),range(F_Best(IDX(ii))));
fprintf('Alpha:  mean = %6.1f,  median = %6.1f, range = %6.1f \n', mean(W(IDX(ii))), median(W(IDX(ii))), range(W(IDX(ii))));
fprintf('Beta:   mean = %6.1f,  median = %6.1f, range = %6.1f \n', mean(C1C2(IDX(ii))), median(C1C2(IDX(ii))), range(C1C2(IDX(ii))));
fprintf('Pop:    mean = %6.1f,  median = %6.1f, range = %6.1f \n', mean(Part(IDX(ii))), median(Part(IDX(ii))), range(Part(IDX(ii))));
