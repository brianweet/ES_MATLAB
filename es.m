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
    mu = 5;
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
    
    % Initialize fitness
    for i = 1:mu
        fp(i,1) = feval(fitnessfct,xp(i,:));
    end
    
    evalcount = evalcount + mu;
    
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

            %TODO BwE: create setting for recombination  method
            
            %intermediate with two parents   
            r1 = select_random(xp);
            r2 = select_random(xp);  
            offspring(i,:) = recombine_intermediate([r1;r2]);

%             %discrete with two parents
%             r1 = select_random(xp);
%             r2 = select_random(xp);            
%             offspring(i,:) = recombine_discrete([r1;r2]);
%
%             %intermediate with all parents   
%             r1 = select_random(xp);
%             r2 = select_random(xp); 
%             offspring(i,:) = recombine_intermediate(xp);
%             
%             %discrete with all parents          
%             offspring(i,:) = recombine_discrete(xp);
            
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

%Recombination methods
function recombined_offspring = recombine_discrete(population)
    % recombined_offspring = recombine_discrete(population)
    %   
    %   Can be used for discrete and global discrete recombination 
    %   population is the list of the parents that will be recombined
    %   The order of population is arbitrary
    %   The size of population should be > 1
    %
    % Author: B. Weeteling
    
    population_size = size(population,1);
    length_of_individual = size(population,2);
    selected_positions = randi(population_size,length_of_individual,1);
    
    recombined_offspring = zeros(length_of_individual,1);
    
    for position_index = 1:length_of_individual
        selected_parent = population(selected_positions(position_index),:);
        recombined_offspring(position_index) = selected_parent(1,position_index);
    end
end

function recombined_offspring = recombine_intermediate(population)
    % recombined_offspring = recombine_intermediate(population)
    %   
    %   Can be used for intermediate and global intermediate recombination 
    %   population is the list of the parents that will be recombined
    %   The order of population is arbitrary
    %   The size of population should be > 1
    %
    % Author: B. Weeteling
    
    recombined_offspring = mean(population);
end

%Selection methods
function selected_individual = select_random(population)
% selected_individual = select_random(population)
%
%   Select 1 individual of population 'population'
%
% Author: B. Weeteling
length_of_individual = size(population,1);
random_individual_idx = randi(length_of_individual,1,1);
selected_individual = population(random_individual_idx,:);
end

function selected_individual = select_proportional(population, population_fitness)
% selected_individual = select_proportional(parents, fitness)
%
%   Select 1 individual of population 'population' with fitness
%   values 'population_fitness' using proportional selection.
%
% Author: Johannes W. Kruisselbrink
% Last modified: October 21, 2009

	cumsum_f = cumsum(population_fitness);
	r = sum(population_fitness) * rand();
	i = 1;
	while (r >= cumsum_f(i))
		i = i + 1;
	end
	selected_individual = population(i, :);

end

function selected_individual = select_tournament(population, population_fitness, tournament_size)
% selected_individual = select_tournament(population, population_fitness, tournament_size)
%
%   Select 1 individual of population 'population' with fitness
%   values 'population_fitness' using tournament selection.
%   'tournament_size' is the tournament size
%
% Author: B. Weeteling

    selected_individuals_index = randi(size(population,1),1,tournament_size);
    selected_individuals = population(selected_individuals_index,:);
    selected_fitness = population_fitness(selected_individuals_index,:);
    
    [~,best_index] = sort(selected_fitness);
    selected_individual = selected_individuals(best_index(1,1), :);
end

