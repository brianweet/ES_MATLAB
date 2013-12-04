function [xp, fp, stat] = es(fitnessfct, n, lb, ub, stopeval)
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
    mu = 2;
    % Amount of child individuals
    lambda = 15;
  
    % Initialize population
    xp = initialize_population(n, mu, lb(1,1), ub(1,1));
    % Create variable for fitness of the population
    fp = zeros(mu,1);
    % Initialize one sigma randomly between upper and lowerbound / 6
    % TODO BwE:  check the sigma initialization
    sigma = (lb(1,1) + (ub(1,1)-lb(1,1))* rand(1,1)) / 6;
    
    % Initialize Tau for 1 sigma strategy
    tau = 1/sqrt(n);
    
    % Create evalcount variable
    evalcount = 0;
  
    % Statistics administration
    stat.name = ['(mu,lambda)-ES mu:' mu ' lambda: ' lambda] ;
    stat.evalcount = 0;
    stat.histsigma = zeros(1, stopeval);
    stat.histf = zeros(1, stopeval);
    
    % Evolution cycle
    while evalcount < stopeval
        % Initialize memory for offspring 
        offspring_fitness = zeros(lambda,1);
        offspring = zeros(lambda,n);
        
        % Steps for mu lambda strategy:
        % 1.Recombine    
        % 2.Mutate        
        % 3.Evaluate        
        % 4.Select
        
        for i = 1:lambda
            % Step 1: Recombination 
            [r1, r2] = selection(xp,fp);
            offspring(i) = recombine_discrete(r1,r2);
            % Step 2: Mutation
            offspring(i) = mutate(offspring(i));
        end 
            
%         %BwE Wrong::Generate offspring from parent xp
%         for i = 1:mu
%             % Mutation using sigma
%             xo(i*lambda-1:i*lambda,:) = ...
%                 repmat(xp(mu,:),[lambda 1]) + sigma * randn(lambda, n);
%         end
        
        % Step 3: Evaluate
        % Evaluate offspring using fitnessfct
        for j = 1:lambda
            offspring_fitness(j,1) = feval(fitnessfct,offspring(j,:));
        end
        
        % Increment eval counter
        evalcount = evalcount + lambda;
        
        % Step 4: Select
        % Pick size(mu) best individuals from offspring and set as mu(t+1)
        [sorted_fitness,idx] = sort(offspring_fitness);
        xp = offspring(idx(1:mu),:);
        fp = sorted_fitness(mu);
        
        % Update sigma
        sigma = sigma * exp(tau*rand);
        
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

function recombined_offspring = recombine_discrete(first_parent,second_parent)
    % recombined_individual = recombine_discrete(first_parent,second_parent)
    %
    %   Recombine two parents to create one offspring
    %   first_parent is one of the parents that will be recombined
    %   second_parent is one of the parents that will be recombined
    %   The order of first_parent and second_parent is arbitrary
    %
    % Author: B. Weeteling
    
    % Create random string (prob 0.5 for either 0 or 1). E.g. 111000
    rnd1 = rand(n,1) > 0.5;        
    % Select bits from first parent where rnd1 == 1. E.g. 101110 + 111000 -> 101000
    first_parent = and(first_parent,rnd1);
    % Select bits from first parent where rnd1 == 0  E.g. 111111 + 000111 -> 000111
    second_parent = and(second_parent,~rnd1);
    % Do the actual recombination E.g. 101000 + 000111 -> 101111
    recombined_offspring = first_parent|second_parent;
end

