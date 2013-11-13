function [xp, fp, stat] = es(fitnessfct, n, lb, ub, stopeval)
% [xp, fp, stat] = es(fitnessfct, n, lb, ub, stopeval)
  % Strategy parameters
  successRate = 0;
  c = 0.817; % 0.817 <= c <= 1
  
  % Initialize
  xp = rand(1,n);
  fp = feval(fitnessfct,xp);
  %TODO BwE:  check the sigma initialization
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
    % generate offspring from parent xp 
    xo = xp + sigma * randn(1, n);
    % evaluate xo using fitnessfct
    fo = feval(fitnessfct,xo); 
    evalcount = evalcount + 1;

    % select best and update success-rate
    if(fo < fp)
        successRate = successRate + 1;
        fp = fo;
        xp = xo;
    end
    
    %update stepsize: change sigma every 5th execution (if needed)
    if(mod(evalcount,5) == 0)
        if(successRate > 2)
            sigma = sigma / c;
        else if(successRate < 2)
                sigma = sigma * c;
            end
        end
        successRate = 0;
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

