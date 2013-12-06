function [xp, fp, stat] = es2(fitnessfct, n, lb, ub, stopeval, mu, lambda)
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
    stat.name = ['(' num2str(mu) ',' num2str(lambda) ')-ES with tau & tau\prime '];
    stat.evalcount = 0;
    stat.histsigma = zeros(1, stopeval);
    stat.histf = zeros(1, stopeval);
    
    % Initialize fitness
    for i = 1:mu
        fp(i) = feval(fitnessfct,xp(i,:));
    end
    stat.histf(1:mu)=min(fp);
    evalcount = evalcount + mu;
    
    % Evolution cycle
    while evalcount < stopeval
        tau_prime = 1 / sqrt(2*n);
        
        % Initialize memory for offspring 
        offspring_fitness = zeros(lambda,1);
        offspring_sigma = repmat(mean(sigma), lambda,1);        
        offspring = repmat(mean(xp), lambda,1);
        
        for i = 1:lambda
            tau = 1 / sqrt(2 * sqrt(n));
            offspring_sigma(i,:) = offspring_sigma(i,:) * exp(tau_prime*randn+tau*randn);
            offspring(i,:) = offspring(i,:) + offspring_sigma(i,:) * randn(1, n);
            offspring_fitness(i,1) = feval(fitnessfct,offspring(i,:));
        end
        
        % Steps for mu lambda strategy:
        % 1.Recombine    
        % 2.Mutate        
        % 3.Evaluate        
        % 4.Select
        
%         for i = 1:lambda
%             % Step 1: Recombination 
%             %TODO BwE: create setting for recombination  method
%             
%             %intermediate with all parents 
%             %[offspring(i,:), offspring_sigma(i,:)] = recombine_intermediate(xp, sigma);
%             
%             % Select parent1
%             [p1, sigma_r1] = select_tournament(xp, sigma, fp, 2);
%             
%             % Remove r1 from xp (no similar parents)
%             index = true(1, size(xp, 1));
%             [~,idx_remove] = ismember(p1, xp,'rows');
%             index(idx_remove) = false;
%             
%             % Select parent2
%             [p2, sigma_r2] = select_tournament(xp(index, :), sigma(index, :), fp(index,:), 2); 
%              
%             % Intermediate with two parents 
%             [offspring(i,:), offspring_sigma(i,:)] = recombine_intermediate([p1;p2], [sigma_r1;sigma_r2]);
% % 
% %             %discrete with two parents           
% %             [offspring(i,:), offspring_sigma(i,:)] = recombine_discrete([p1;p2], [sigma_r1;sigma_r2]);
% %           
% %             %discrete with all parents          
% %             [offspring(i,:), offspring_sigma(i,:)] = recombine_discrete(xp, sigma);
%             
%             % Step 2: Mutation
%             [offspring(i,:), offspring_sigma(i,:)] = mutate(offspring(i,:),offspring_sigma(i,:), tau_prime);
%         end 
        
        % Step 3: Evaluate
        % Evaluate offspring using fitnessfct
%         for j = 1:lambda
%             offspring_fitness(j,1) = feval(fitnessfct,offspring(j,:));
%         end
        
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

%Mutation methods
function [mutated_individual, mutated_sigma] = mutate(individual, sigma, tau_prime)
    %   [mutated_offspring, mutated_sigma] = mutate(individual, sigma)
    %   
    %   Can be used to mutate the individual and the sigma 
    %   individual is the individual that will be mutated
    %   sigma is the sigma that will be mutated
    %
    % Author: B. Weeteling
    
    n = length(individual);
    tau = 1 / sqrt(2 * sqrt(n));
    
    
%     mutated_sigma = sigma * exp(tau_prime*randn);
%     mutated_individual =  individual + mutated_sigma * randn(1,n); 
    
    % TODO BwE: understand how to implement this 
    mutated_sigma = sigma * exp(tau_prime * randn + tau * randn);
    mutated_individual = individual + (mutated_sigma .* randn(n, 1))'; 
end

%Recombination methods
function [recombined_offspring, recombined_sigma] = recombine_discrete(population, sigma)
    %  [recombined_offspring, recombined_sigma]  = recombine_discrete(population, sigma)
    %   
    %   Can be used for discrete and global discrete recombination 
    %   population is the list of the parents that will be recombined
    %   sigma is the list of sigma values for the parent population
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
    
    recombined_sigma = mean(sigma);
end

function [recombined_offspring, recombined_sigma] = recombine_intermediate(population, sigma)
    % recombined_offspring = recombine_intermediate(population)
    %   
    %   Can be used for intermediate and global intermediate recombination 
    %   population is the list of the parents that will be recombined
    %   sigma is the list of sigma values for the parent population
    %   The order of population is arbitrary
    %   The size of population should be > 1
    %
    % Author: B. Weeteling
    
    recombined_offspring = mean(population);
    recombined_sigma = mean(sigma);
end

%Selection methods
function [selected_individual, selected_sigma] = select_random(population, sigma)
% selected_individual = select_random(population)
%
%   Select 1 individual of population 'population'
%
% Author: B. Weeteling
length_of_individual = size(population,1);
random_individual_idx = randi(length_of_individual,1,1);
selected_individual = population(random_individual_idx,:);
selected_sigma = sigma(random_individual_idx,:);
end

function [selected_individual, selected_sigma] = select_proportional(population, sigma, population_fitness)
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
    selected_sigma = sigma(i,:);

end

function [selected_individual, selected_sigma] = select_tournament(population, sigma, population_fitness, tournament_size)
% selected_individual = select_tournament(population, population_fitness, tournament_size)
%
%   Select 1 individual of population 'population' with fitness
%   values 'population_fitness' using tournament selection.
%   'tournament_size' is the tournament size
%
% Author: B. Weeteling
    if tournament_size > size(population,1) 
        error('tournament size should be smaller then population size'); 
    end
    
    selected_individuals_index = randperm(size(population,1),tournament_size);
    selected_individuals = population(selected_individuals_index,:);
    selected_fitness = population_fitness(selected_individuals_index,:);
    selected_sigmas = sigma(selected_individuals_index,:);
    
    [~,best_index] = sort(selected_fitness);
    selected_individual = selected_individuals(best_index(1,1), :);
    selected_sigma = selected_sigmas(best_index(1,1), :);
end

