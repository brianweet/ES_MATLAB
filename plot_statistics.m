function [] = plot_statistics(stat, fitnessfct, n, lb, ub, stopeval, runs)
    endresult = 0;  
    for i = 1 : runs
        histf(:,i) = stat(i).histf(1:stopeval);
        endresult = endresult+ stat(i).histf(stopeval);
    end
    endresult = endresult / runs;
    disp(num2str(endresult))
    figure;
    tmp = mean(histf, 2);
    tmp(stopeval)
    plot(mean(histf, 2))
    title(stat(1).name);
    xlabel('evaluations')
    ylabel('fitness')
end

