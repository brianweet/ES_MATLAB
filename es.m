function [xp, fp, stat] = es(fitnessfct, n, lb, ub, stopeval)
  
    % Strategy parameters
    % Amount of parent individuals
    mu = 2;
    % Amount of child individuals
    lambda = 2;
  
    % Initialize population
    xp = initialize_population(n, mu, lb(1,1), ub(1,1));
    % Create variable for fitness of the population
    fp = zeros(mu,1);
    % Initialize one sigma randomly between upper and lowerbound / 6
    % TODO BwE:  check the sigma initialization
    sigma = (lb(1,1) + (ub(1,1)-lb(1,1))* rand(1,1)) / 6;
    % Create evalcount variable
    evalcount = 0;
    % Set empty succes rate array (we're keeping track of every evaluation,not needed per se)
    successRate = zeros(stopeval,1);
  
    % Statistics administration
    stat.name = ['(mu+lambda)-ES mu:' mu ' lambda: ' lambda] ;
    stat.evalcount = 0;
    stat.histsigma = zeros(1, stopeval);
    stat.histf = zeros(1, stopeval);
    
    % Evolution cycle
    while evalcount < stopeval
        % Initialize offspring fitness var
        fo = zeros(mu*lambda,1);
        xo = zeros(mu*lambda,n);
        
        % Generate offspring from parent xp
        for i = 1:mu
            xo(i*lambda-1:i*lambda,:) = ...
                repmat(xp(mu,:),[lambda 1]) + sigma * randn(lambda, n);
        end
        
        % Evaluate new population using fitnessfct
        for j = 1:mu*lambda
            fo(j,1) = feval(fitnessfct,xo(j,:));
        end
        
        xp = xo;
        fp = fo;
        
        % Increment eval counter
        evalcount = evalcount + mu*lambda;
        
        % Update stepsize: change sigma every n-th execution (if needed)
%         if mod(evalcount-1,n) == 0
%             startIndex = 1;
%             if evalcount > (10*n)
%                 startIndex = evalcount-10*n;
%             end
%             %Ps is the relative frequency of successful mutations measured over
%             % 10 * n trials (slides 5.3)
%             Ps = sum(successRate(startIndex:evalcount,1))/(evalcount-startIndex);
%             if(Ps > .2)
%                 sigma = sigma / c;
%             else if(Ps < .2)
%                     sigma = sigma * c;
%                 end
%             end
%         end
        
        % Statistics administration
        stat.histsigma(evalcount) = sigma;% stepsize history
        stat.histf(evalcount) = min(fp);% fitness history

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

function xp = initialize_population(n, mu, lower, upper)
    % xp = initialize_population(n, mu)
    %
    %   Initialize new population with values within a given range
    %   n is the number of dimensions
    %   mu is the number of individuals
    %   lower is the lowerbound
    %   upper is the upperbound
    %
    % Author: B. Weeteling
    xp = lower+(upper-lower).*rand(mu,n);
end

