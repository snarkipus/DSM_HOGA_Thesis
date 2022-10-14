function [J,G_init] = HP_init(P)

%% nth-Order High-Pass NTF Objective Function Initialization
%  //////////////////////////////////////////////////////////////////////////
% ///
% ///       Usage: [F_s,G_init] = HP_init(P)
% ///
% ///   Arguments: [int]   P: Population Size (# of chromosomes)
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
% /// File Name: HP_init.m
% ///
% /// Description:
% /// M-File which initializes the HTGA for realizing a high-pass NTF.
% ///
% /// Author: M. Jackson
% ///
% /// Changelog:
% ///
% ///    [DATE]      [VERSION]
% ///
% ///    09/04/07    0.1
% ///    generalized n'th order objective function initialization adapted from
% ///    HP4_init.m (static 4th order function)
% ///
% ///    09/05/07    0.11
% ///    -  added regional seed constraints
% ///    -  changed seeding to be sane in terms of coefficient dependence
% ///
% ///    09/05/07    0.2 (MAJOR REWORK)
% ///    -  changed seeding to chebyshev polynomial with normally dithered
% ///       coefficients 
% ///
% ///    09/08/07    0.21
% ///    -  changed cheby2() cut-off frequency to reflect discrete time filter
% ///       normalization wrt 1 vs. pi
% ///    -  made 'order' a global variable
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
theta  = zeros(1,2*M+4*N+1);
G_init = zeros(2*M+4*N+1,P);

%% Seeding Constraints

% Design a 'seed' chebychev II filter
[seed.num,seed.den] = cheby2(order,DR,1/OSR,'high');

% Convert the 'seed' filter to a product of second-order sections (SOS)
[seed.sos,seed.g] = tf2sos(seed.num,seed.den); 

% Translate the 'seed' SOS into a parametric chromosome (theta)
i = 1;
j_start = 1;
j_stop  = 0;

if( M==1 ) % odd-order requires a first order section
    
    seed.theta(1) = seed.sos(1,2);
    seed.theta(2) = seed.sos(1,5);
    i = 3;
    j_start = 2;
    j_stop  = 1;
    
end
    
% 2nd-Order Numerator Coefficients
for j = j_start : (N + j_stop)
    
    seed.theta(i)   = seed.sos(j,2);
    seed.theta(i+1) = seed.sos(j,3);
    i = i +2;
    
end

% 2nd-Order Denominator Coefficients
for j = j_start : (N + j_stop)
    
    seed.theta(i)   = seed.sos(j,5);
    seed.theta(i+1) = seed.sos(j,6);
    i = i +2;
    
end

% Gain
seed.theta(2*M + 4*N + 1) = seed.g;

% Generate a random population
G_init = 1e-2*randn(2*M+4*N+1,P);

% Seed the target filter;
for i = 1:P
    G_init(1:end,i) = G_init(1:end,i) + seed.theta';
end

% Objective Function Handle
J = @HP_obj;

end
