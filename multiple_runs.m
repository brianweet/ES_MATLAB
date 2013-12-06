function [] = multiple_runs(FUN, fitnessfct, n, lb, ub, stopeval, runs, mu, lambda)
  for i = 1 : runs
    [xopt, fopt, stat(i)] = feval(FUN,fitnessfct, n, lb, ub, stopeval, mu, lambda);
  end
  save('statistics.mat', 'stat')
end

