% es(@sphere, 10, -10, 10, 10000)

% multiple_runs(@sphere, 10, -10, 10, 10000, 10)
% load('statistics.mat');
% plot_statistics(stat, @sphere, 10, -10, 10, 10000, 10)

multiple_runs(@ackley, 10, -10, 10, 10000, 10)
load('statistics.mat');
plot_statistics(stat, @ackley, 10, -10, 10, 10000, 10)
 
% multiple_runs(@rosenbrock, 10, -10, 10, 10000, 10)
% load('statistics.mat');
% plot_statistics(stat, @rosenbrock, 10, -10, 10, 10000, 10)

% [xopt, fopt] = lastname1_lastname2_es(@bbf1, 30, -100 * ones(1,30), 
% 100 * ones(1,30), 10000) 

% multiple_runs(@bbf1, 30, -100 * ones(1,30), 100 * ones(1,30), 10000,10)
% load('statistics.mat');
% plot_statistics(stat, @bbf1, 30, -100 * ones(1,30), 100 * ones(1,30), 10000, 10)
