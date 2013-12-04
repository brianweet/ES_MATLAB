% es(@sphere, 10, -10, 10, 10000)

% multiple_runs(@es,@sphere, 10, -10, 10, 10000, 10)
% load('statistics.mat');
% plot_statistics(stat, @sphere, 10, -10, 10, 10000, 10)

% multiple_runs(@es,@ackley, 10, -10, 10, 10000, 10)
% load('statistics.mat');
% plot_statistics(stat, @ackley, 10, -10, 10, 10000, 10)
 
% multiple_runs(@es,@rosenbrock, 10, -10, 10, 10000, 10)
% load('statistics.mat');
% plot_statistics(stat, @rosenbrock, 10, -10, 10, 10000, 10)

% [xopt, fopt] = lastname1_lastname2_es(@bbf1, 30, -100 * ones(1,30), 
% 100 * ones(1,30), 10000) 

multiple_runs(@es, @bbf1, 30, -100 * ones(1,30), 100 * ones(1,30), 10000, 1)
%load('statistics.mat');
%plot_statistics(stat, @bbf1, 30, -100 * ones(1,30), 100 * ones(1,30), 10000, 1)

% multiple_runs(@es,@bbf2, 30, -100 * ones(1,30), 100 * ones(1,30), 10000,10)
% load('statistics.mat');
% plot_statistics(stat, @bbf2, 30, -100 * ones(1,30), 100 * ones(1,30), 10000, 10)
% 
% 
% multiple_runs(@es,@bbf3, 30, -100 * ones(1,30), 100 * ones(1,30), 10000,10)
% load('statistics.mat');
% plot_statistics(stat, @bbf3, 30, -100 * ones(1,30), 100 * ones(1,30), 10000, 10)
% 
% 
% multiple_runs(@es,@bbf4, 30, -100 * ones(1,30), 100 * ones(1,30), 10000,10)
% load('statistics.mat');
% plot_statistics(stat, @bbf4, 30, -100 * ones(1,30), 100 * ones(1,30), 10000, 10)
% 
% 
% multiple_runs(@es,@bbf5, 30, -100 * ones(1,30), 100 * ones(1,30), 10000,10)
% load('statistics.mat');
% plot_statistics(stat, @bbf5, 30, -100 * ones(1,30), 100 * ones(1,30), 10000, 10)

