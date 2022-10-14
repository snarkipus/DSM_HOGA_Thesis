function [J,G_init] = LPSTF_init(P)

%% nth-Order Low-Pass STF Objective Function Initialization
%  //////////////////////////////////////////////////////////////////////////
% ///
% ///       Usage: [F_s,G_init] = LPSTF_init(P,J)
% ///
% ///   Arguments: [int]   P:              Population Size (# of chromosomes)
% ///              [fcn]   J:              Objective Function Handle
% ///
% /// Returns:     [stuct] F_s:            Fitness Structure
% ///              [mat]   G_init:         Initial Population Matrix
% ///
% /// Definitions: <F_s>   F_s.fitness:    Fitness Vector
% ///                      F_s.best_chr:   Best Chromosome
% ///                      F_s.best_cost:  Best Fitness Value
% ///                      F_s.worst_chr:  Worst Chromosome
% ///                      F_s.worst_cost: Worst Fitness Value
% ///
% /////////////////////////////////////////////////////////////////////////

%% File Information
%  //////////////////////////////////////////////////////////////////////////
% ///
% /// File Name: LPSTF_init.m
% ///
% /// Description:
% /// M-File which initializes the HTGA for realizing a low-pass STF.
% ///
% /// Author: M. Jackson
% ///
% /// Changelog:
% ///
% ///    [DATE]      [VERSION]
% ///
% ///    09/09/07    0.1 (initial work)
% ///
% /////////////////////////////////////////////////////////////////////////

%% Variable Initialization

global OSR DR order

% Process the filter order
if( mod(order,2)==1 )
    
    % odd-order filter requires 1st-order section
    N = (order-1)/2;    % number of 2nd-order sections
    M = 1;              % number of 1st-order sections 
    
else
    
    % all 2nd-order sections
    N = order/2;
    M = 0;
    
end

% Preallocation
gamma  = zeros(1,M+2*N+1);
G_init = zeros(M+2*N+1,P);

%% Seeding Constraints

% Design a 'seed' chebychev II filter
[seed.num,seed.den] = cheby2(order,DR,1/OSR,'low');

% Convert the 'seed' filter to a product of second-order sections (SOS)
[seed.sos,seed.g] = tf2sos(seed.num,seed.den); 

% Translate the 'seed' SOS into a parametric chromosome (theta)
i = 1;

if( M==1 ) % odd-order requires a first order section
    
    seed.gamma(1) = seed.sos(1,2);
    i = 2;
    
end
    
% 2nd-Order Numerator Coefficients
for j = 1:N
    
    seed.gamma(i)   = seed.sos(j,2);
    seed.gamma(i+1) = seed.sos(j,3);
    i = i +2;
    
end

% Gain
seed.gamma(M+2*N+1) = seed.g;

% Generate a random population
G_init = 1e-2*randn(M+2*N+1,P);

% Seed the target filter;
for i = 1:P
    G_init(1:end,i) = G_init(1:end,i) + seed.gamma';
end

%% Initial Population Fitness Evaluation
J = @LPSTF_obj;

end
