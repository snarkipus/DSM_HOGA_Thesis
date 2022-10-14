function G_pool = LR_Selection(G, J)

%% Linear-Ranking Selection
%  //////////////////////////////////////////////////////////////////////////
% ///
% ///       Usage: G_pool = LR_Selection(G, J)
% ///
% ///   Arguments: [mat]      G: Generation
% ///              [fcn]      J: Objective Function Handle
% ///
% ///     Returns: [mat] G_pool: Selected Generation Pool
% ///
% /// Definitions: NONE
% ///
% //////////////////////////////////////////////////////////////////////////

%% File Information
%  //////////////////////////////////////////////////////////////////////////
% ///
% /// File Name: LR_Selection.m
% ///
% /// Description:
% /// Linear Ranking Selection [Baker, 1987].
% ///
% /// Author: M. Jackson
% ///
% /// Changelog:
% ///
% ///    [DATE]      [VERSION]
% ///
% ///    08/19/07    0.1 (initial work - HTGAv2.m)
% ///
% ///    09/01/07    0.2
% ///    extracted from HTGAv3.m to create stand-alone LR_Selection.m
% ///
% ///    09/02/07    0.21
% ///    fixed header information
% ///
% //////////////////////////////////////////////////////////////////////////

%% Initialization

% Population Information
P = size(G,2);  % Number of Chromosomes
N = size(G,1);  % Number of Traits per Chromosome

% Preallocation
prob.vec  =   zeros(1,P);   % Selection Probability Vector
prob.dist =   zeros(1,P);   % Selection Probability Distribution
Gsort     = zeros(N+1,P);   % temporary sorting matrix

%% Fitness Evaluation
% *TODO:* Remove this redundant evaluation

% Evaluate Population Fitness
F_s = Eval_Fitness(G, J);

% Process Arguments
F          = F_s.fitness;
worst_cost = F_s.worst_cost;
best_chr   = F_s.best_chr;

%% Linear-Ranking Scheme
% *TODO:* Paramaterize the min/max pdf's to allow full range of selection pressure

% Invert / Shift Fitness Vector (non-negative)
F = (-F) + max(F_s.fitness); 

% Sort Population Based on Fitness
Gsort  = [F' G']';              % Append fitness values to population matrix
Gsort  = sortrows(Gsort',1)';   % Sort population
G      = Gsort;                 % Reassign population
G(1,:) = [];                    % Delete fitness row

 % Pre-select the best progenitor (Elitism)
G_pool(:,1) = best_chr;

% Calculate the Selection Probability for each chromosome (linearly ranked)
eta_min = 0.9;
eta_max = 2 - eta_min;

for i = 1:P
    prob.vec(i) = (1/P)*(eta_min + (eta_max-eta_min)*((i-1)/(P-1)));
end

prob.dist = cumsum(prob.vec); % probability distribution vector

% Select the survivors from the previous population
new_idx = 2;

for i = 1:P

    fit_rand = rand();

    if( prob.dist(i) > fit_rand )
        G_pool(:,new_idx) = G(:,i);
        new_idx = new_idx + 1;
    end

end

end