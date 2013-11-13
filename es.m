function [xp, fp, stat] = es(fitnessfct, n, lb, ub, stopeval)
% [xp, fp, stat] = es(fitnessfct, n, lb, ub, stopeval)
  % Strategy parameters
  ...
  tau = 1/sqrt(n);
  % Initialize
  xp = rand(1,n);
  fp = feval(fitnessfct,xp);
  sigma = (lb + (ub-lb).* rand(1,1) ) / 6;
  evalcount = 0;

  % Statistics administration
  stat.name = '(1+1)-ES';
  stat.evalcount = 0;
  stat.histsigma = zeros(1, stopeval);
  stat.histf = zeros(1, stopeval);

  % Evolution cycle
  while evalcount < stopeval
    % Generate offspring and evaluate
    sigmao = sigma * exp(tau*randn);
    xo = xp + sigmao * randn(1, n); % generate offspring from parent xp 
    fo = feval(fitnessfct,xo); % evaluate xo using fitnessfct
    evalcount = evalcount + 1;

    % select best and update success-rate and update stepsize
    % Important: MINIMIZATION!
    if(fo < fp)
        fp = fo;
        xp = xo;
        sigma = sigmao;
    end
    
    % Statistics administration
    stat.histsigma(evalcount) = sigma;% stepsize history
    stat.histf(evalcount) = fp;% fitness history

%     % if desired: plot the statistics
%     % Plot statistics
%     clf
%     subplot(2,1,1)
%     plot(stat.histf(1:evalcount))
%     subplot(2,1,2)
%     plot(stat.histsigma(1:evalcount))
%     drawnow()
  end

end

