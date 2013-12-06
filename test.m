function ParentPop = test()
mu = 5;                 % number of parents
lambda = 15;            % number of offspring
yInit = ones(30,1);     % initial parent vector 
sigmaInit = 1;          % initial global mutation strength sigma 
sigmaMin = 1e-10;       % ES stops when sigma is smaller than sigmaMin

% initialization:
n = length(yInit);      % determine search space dimensionality n   
tau = 1/sqrt(2*n);      % self-adaptation learning rate
% initializing individual population:
Individual.y = yInit;
Individual.sigma = sigmaInit;
Individual.F = fitness(Individual.y);
for i=1:mu; ParentPop{i} = Individual; end


count = 0;
while(count < 10000)
 Recombinant = recombine(ParentPop);              % recombine parents
 for l = 1:lambda;                                % generate lambda offspring
  OffspringIndividual.sigma = Recombinant.sigma * exp(tau*randn); % mutate sigma
  OffspringIndividual.y = Recombinant.y + OffspringIndividual.sigma * randn(n, 1); % mutate object parameter
  OffspringIndividual.F = fitness(OffspringIndividual.y); % determine fitness
  OffspringPop{l} = OffspringIndividual;                  % offspring complete
 end;
 ParentPop = SortPop(OffspringPop, mu);   % sort population
 disp(ParentPop{1}.F);                    % display best fitness in population
 if ( ParentPop{1}.sigma < sigmaMin ) 
     break 
 end % termination criterion
 count = count+lambda;
end
end


% function to be optimized (sphere test function as an example, to be changed):
function out = fitness(x)
    x=x';
    out = feval(@bbf5,x);
    %out = sum(x.*x); 
end

% this sorts the population according to the individuals' fitnesses:
function sorted_pop = SortPop(pop, mu);
 for i=1:length(pop); fitnesses(i) = pop{i}.F; end;
 [sorted_fitnesses, index] = sort(fitnesses);
 for i=1:mu; sorted_pop{i} = pop{index(i)}; end
end

% this performs intermediate (multi-) recombination: 
function r = recombine(pop)
 r.sigma = 0;
 r.y = 0;
 for i=1:length(pop) 
  r.sigma = r.sigma + pop{i}.sigma;
  r.y = r.y + pop{i}.y;
 end;
 r.sigma = r.sigma/length(pop); 
 r.y = r.y/length(pop);
end

