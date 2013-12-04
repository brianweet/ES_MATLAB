function [] = multiple_runs(FUN, fitnessfct, n, lb, ub, stopeval, runs)
  for i = 1 : runs
    [xopt, fopt, stat(i)] = feval(FUN,fitnessfct, n, lb, ub, stopeval);
  end
  save('statistics.mat', 'stat')
end

