function legendname = plot_statistics(stat, fitnessfct, n, lb, ub, stopeval, runs, linestyle)
    
    legendname = stat(1).name;
    valuesAt1000 = zeros(runs,1);
    endvalues = zeros(runs,1);
    
    for i = 1 : runs
        histf(:,i) = stat(i).histf(1:stopeval);
        
        endvalues(i) = stat(i).histf(10000);
        valuesAt1000(i) = stat(i).histf(1000);
    end
    
    disp(legendname)
    disp(['1000 ' ...
          num2str(min(valuesAt1000)) ' ' ... 
          num2str(max(valuesAt1000)) ' ' ...
          num2str(mean(valuesAt1000)) ' ' ...
          num2str(std(valuesAt1000))])
        
    disp(['10000 ' ...
          num2str(min(endvalues)) ' ' ... 
          num2str(max(endvalues)) ' ' ...
          num2str(mean(endvalues)) ' ' ...
          num2str(std(endvalues))])
    disp(' ')
    
    tmp = mean(histf, 2);
    tmp(stopeval);
    plot(mean(histf, 2),linestyle)
    xlabel('evaluations')
    ylabel('fitness')
end

