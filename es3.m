function [xp, fp, stat] = es3(fitnessfct, n, lb, ub, stopeval, mu, lambda)
% [xp, fp, stat] = es(fitnessfct, n, lb, ub, stopeval)
%
%   Run mu+lambda ES with the given objective fitness function
%   n is the number of dimensions
%   lb is the lowerbound of the search space
%   ub is the upperbound of the search space
%   stopeval is the max amount of function evaluations
%
% Author: B. Weeteling
  
    % Strategy parameters
    % Amount of parent individuals
    %mu = 5;
    % Amount of child individuals
    %lambda = 12;
  
    % Initialize population
    xp = initialize_population(n, mu, lb(1,1), ub(1,1));
    % Create variable for fitness of the population
    fp = zeros(mu,1);
    
    % Initialize one sigma randomly between upper and lowerbound / 6
    % TODO BwE:  check the sigma initialization
    sigma = (lb(1,1) + (ub(1,1)-lb(1,1)) .* rand(mu,1))/6;
    
    % Create evalcount variable
    evalcount = 0;
  
    % Statistics administration
    stat.name = ['(' num2str(mu) ',' num2str(lambda) ')-ES only tau\prime'];
    stat.evalcount = 0;
    stat.histsigma = zeros(1, stopeval);
    stat.histf = zeros(1, stopeval);
    
    % Initialize fitness
    for i = 1:mu
        fp(i) = feval(fitnessfct,xp(i,:));
    end
    stat.histf(1:mu)=fp;
    evalcount = evalcount + mu;
    
    % Evolution cycle
    while evalcount < stopeval
        % Steps for mu lambda strategy:
        % 1.Recombine   
        % 2.Mutate        
        % 3.Evaluate        
        % 4.Select
        
        tau_prime = 1 / sqrt(2*n);
        
        % Initialize memory for offspring 
        offspring_fitness = zeros(lambda,1);
        % 1.Recombine
        offspring_sigma = repmat(mean(sigma), lambda,1);        
        offspring = repmat(mean(xp), lambda,1);
        
        for i = 1:lambda
            % 2.Mutate
            offspring_sigma(i,:) = offspring_sigma(i,:) * exp(tau_prime*randn);
            offspring(i,:) = offspring(i,:) + offspring_sigma(i,:) * randn(1, n);
            % 3.Evaluate
            offspring_fitness(i,1) = feval(fitnessfct,offspring(i,:));
        end
        
        % Increment eval counter
        evalcount = evalcount + lambda;
        
        % Step 4: Select
        % Pick size(mu) best individuals from offspring and set as mu(t+1)
        [sorted_fitness,idx] = sort(offspring_fitness);
        xp = offspring(idx(1:mu),:);
        fp = sorted_fitness(1:mu,:);
        sigma = offspring_sigma(idx(1:mu),:);
        
        % Statistics administration
        stat.histf(evalcount-lambda:evalcount) = min(fp);    % fitness history
    end
end

function population = initialize_population(n, mu, lowerbound, upperbound)
    % xp = initialize_population(n, mu)
    %
    %   Initialize new population with values within a given range
    %   n is the number of dimensions
    %   mu is the number of individuals
    %   lowerbound is the lowerbound
    %   upperbound is the upperbound
    %
    % Author: B. Weeteling
    population = lowerbound+(upperbound-lowerbound).*rand(mu,n);
end

