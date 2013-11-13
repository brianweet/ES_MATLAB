%es(@sphere, 10, -10, 10, 10000)

multiple_runs(@sphere, 10, -10, 10, 10000, 20)
load('statistics.mat');
plot_statistics(stat, @sphere, 10, -10, 10, 10000, 20)