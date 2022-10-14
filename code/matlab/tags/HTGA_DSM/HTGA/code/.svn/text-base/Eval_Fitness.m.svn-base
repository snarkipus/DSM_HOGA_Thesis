function F_s = Eval_Fitness(G, J)

%//////////////////////////////////////////////////////////////////////////
%///
%/// F_s = Eval_Fitness(G, J)
%///
%/// Description:
%/// Evaluate the fitness of the population
%///
%/// Arguments:     [mat]   G:              Population Matrix
%///                [fcn]   J:              Fitness Function Handle
%///
%/// Returns:       [stuct] F_s:            Fitness Structure
%///
%/// Definitions:   <F_s>   F_s.fitness:    Fitness Vector
%///                        F_s.best_chr:   Best Chromosome
%///                        F_s.best_cost:  Best Fitness Value
%///                        F_s.worst_chr:  Worst Chromosome
%///                        F_s.worst_cost: Worst Fitness Value
%///
%//////////////////////////////////////////////////////////////////////////

% Population Information
P = size(G,2);  % Number of Chromosomes
N = size(G,1);  % Number of Traits per Chromosome

% Initialization
best_chr   = zeros(1,N);
best_cost  =    realmax;
worst_chr  = zeros(1,N);
worst_cost =   -realmax;
F          = zeros(1,P);

for i = 1:P

    F(i) = J(G(:,i));

    % Retain the 'Best' Information
    if(F(i) < best_cost)
        best_chr  = G(:,i);
        best_cost = F(i);
    end

    % Retain the 'Worst' Information
    if(F(i) > worst_cost)
        worst_chr  = G(:,i);
        worst_cost = F(i);
    end

end

F_s.fitness    = F;
F_s.best_chr   = best_chr;
F_s.best_cost  = best_cost;
F_s.worst_chr  = worst_chr;
F_s.worst_cost = worst_cost;

end
