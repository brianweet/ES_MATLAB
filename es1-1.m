function [xp, fp, stat] = es(fitnessfct, n, lb, ub, stopeval)
% [xp, fp, stat] = es(fitnessfct, n, lb, ub, stopeval)
  % Strategy parameters
  successRate = zeros(stopeval,1);
  c = 0.817; % 0.817 <= c <= 1
  
  % Initialize
  xp = rand(1,n);
  fp = feval(fitnessfct,xp);
  %TODO BwE:  check the sigma initialization
  sigma = (lb(1,1) + (ub(1,1)-lb(1,1))* rand(1,1) ) / 6;
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
        successRate(evalcount,1) = 1;
        fp = fo;
        xp = xo;
    end
    
    %update stepsize: change sigma every n-th execution (if needed)
    if mod(evalcount-1,n) == 0
        startIndex = 1;
        if evalcount > (10*n)
            startIndex = evalcount-10*n;
        end
        %Ps is the relative frequency of successful mutations measured over
        % 10 * n trials (slides 5.3)
        Ps = sum(successRate(startIndex:evalcount,1))/(evalcount-startIndex);
        if(Ps > .2)
            sigma = sigma / c;
        else if(Ps < .2)
                sigma = sigma * c;
            end
        end
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

