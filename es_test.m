handles = {@bbf1, @bbf2, @bbf3, @bbf4, @bbf5};
runs = 20;

mu = 3;
lambda = 12;
for i = 1:length(handles)
    handle = handles{i};
    func2str(handle);
    
    figure;
    hold on
    multiple_runs(@es2, handle, 30, -100 * ones(1,30), 100 * ones(1,30), 10000, runs, mu, lambda);
    load('statistics.mat');
    legendname1 = plot_statistics(stat, handle, 30, -100 * ones(1,30), 100 * ones(1,30), 1000, runs, 'r');
    
    multiple_runs(@es2, handle, 30, -100 * ones(1,30), 100 * ones(1,30), 10000, runs, 15, 100);
    load('statistics.mat');
    legendname2 = plot_statistics(stat, handle, 30, -100 * ones(1,30), 100 * ones(1,30), 1000, runs, 'g');

    multiple_runs(@es1_1, handle, 30, -100 * ones(1,30), 100 * ones(1,30), 10000, runs, mu, lambda);
    load('statistics.mat');
    legendname3 = plot_statistics(stat, handle, 30, -100 * ones(1,30), 100 * ones(1,30), 1000, runs, ':b');

    legend(legendname1, legendname2, legendname3)
    hold off
    
end