function Darwin_s = HOGA(config_s,test_s)

%% Hybrid Orthogonal Genetic Algorithm (HOGA) 
%  //////////////////////////////////////////////////////////////////////////
% ///
% ///       Usage: Darwin_s = HOGA(J)
% ///
% ///   Arguments: [struct] test_s   : Objective Function
% ///              [struct] config_s : HOGA configuration
% ///
% ///     Returns: [struct] Darwin_s : Evolved Results 
% ///
% /// Definitions:
% ///
% ///    <Darwin_s>      [struct]    Darwin.F_s:          Fitness Structure
% ///     |              [struct]    Darwin.metric_s:     Metrics Structure
% ///     |
% ///     |-------<F_s>  [vec]       F_s.best_chr:       Best Chromosome
% ///     |              [dbl]       F_s.best_cost:      Best Fitness Value
% ///     |              [vec]       F_s.worst_chr:      Worst Chromosome
% ///     |              [dbl]       F_s.worst_cost:     Worst Fitness Value
% ///     |              [vec]       F_s.fitness:        Fitness Vector
% ///     | 
% ///     |- <metric_s>  [vec]       metric_s.avg_cost:  Average Cost 
% ///                    [vec]       metric_s.var_track: Fitness Variance
% ///                    [vec]       metric_s.best_cost: Best Fitness
% ///                    [vec]       metric_s.worst_cost Worst Fitness   
% ///
% /// File Dependencies:
% ///  - Crossover.m:      Single-Point Crossover Operator
% ///  - Eval_Fitness.m:   Fitness Evaluation Routine
% ///  - LR_Selection.m:   Linear-Ranking Selection Operator
% ///  - Mutate.m:         Point-Source Mutation Operator
% ///  - taguchi_method.m: Taguchi Method Crossover Operator
% ///  - Setup_Obj_Fcn.m:  Objective Function Library File
% ///
% //////////////////////////////////////////////////////////////////////////

%% File Information
%  //////////////////////////////////////////////////////////////////////////
% ///
% /// File Name: HTGA.m
% ///
% /// Description:
% /// Script for prototyping and debugging of the HTGA [Tsai et al,2007].
% ///
% /// Author: M. Jackson
% ///
% /// Changelog:
% ///
% ///    [DATE]      [VERSION]
% ///
% ///    08/07/07    0.1 (initial work)
% ///
% ///    08/17/07    0.11
% ///    implimented initialization and example objective function
% ///
% ///    08/18/07    0.12
% ///    implemented roulette-wheel based selection [DeJong,1975]
% ///
% ///    08/19/07    0.13
% ///    changed roulette-wheel based selection to stochastic universe sampling
% ///    [Baker, 1989]
% ///
% ///    08/21/07    0.2
% ///    implemented one-cut point crossover and mutation routines
% ///
% ///    08/22/07    0.21
% ///    restructured init, select, and mutation
% ///
% ///    08/24/07    0.22
% ///    fixed the mutation routine
% ///
% ///    08/31/07    0.3 (MAJOR RE-WORK)
% ///     -  made functions independent of strict sizing
% ///     -  added additional sorting stage to restrict population size to
% ///        only the initial 'P' (truncation selection)
% ///     -  fixed early convergence by removing zeros() preallocation
% ///     -  changed file name from 'HTGAv2' to 'HTGA_DSM'
% ///     -  started developing IIR objective function
% ///
% ///    09/01/07    0.31
% ///     -  created project directory $HOME/work/HTGA
% ///     -  split development script into functional code modules
% ///     -  developed initial 2nd-order HP filter ojective function
% ///     -  changed file name from 'HTGA_NTF' to 'HTGA_NTF'
% ///
% ///    09/02/07    0.32
% ///     -  developed initial 4th-order HP filter objective function
% ///     -  continued to polish individual code modules
% ///
% ///    09/03/07    0.33
% ///    successfully evolved equiripple HP/LP 4th order filters
% ///
% ///    09/04/07    0.4
% ///     -  successfully implemented the Taguchi Method
% ///     -  generalized the HP init and objective functions for n'th order
% ///
% ///    09/08/07    0.41
% ///     -  generalized the HTGA to accept funcion handles to the objective
% ///        and initialization function
% ///     -  changed the filename to 'HTGA.m'
% ///     -  removed DEBUG functions from main body
% ///
% ///    05/11/08    0.5
% ///     -  changed the name from 'HTGA' to 'HOGA'
% ///     -  Generalized the algorithm for ICSENG testing
% ///
% //////////////////////////////////////////////////////////////////////////

%%  GA Parameter Initialization

% tic

global func_evals;

% HOGA Configuration
P  = config_s.P;
CR = config_s.CR; 
MR = config_s.MR;

max_iter     = config_s.max_iter;
min_iter     = config_s.min_iter;
converge_min = config_s.converge_min;
epsilon      = config_s.epsilon;

%% Main( )

%//////////////////////////////////////////////////////////////////////////
%/// INITIALIZATION

converge  = 0;
loopcount = 1;
best_cost = realmax;

[J,G] = Setup_Obj_Fcn(P, test_s.N, test_s.type);
F_s   = Eval_Fitness(G,J);
func_evals = func_evals + P;
    
%//////////////////////////////////////////////////////////////////////////
%/// MAIN LOOP
while( (converge < converge_min + min_iter) )
    
    %//////////////////////////////////////////////////////////////////////
    %/// SELECTION
    G_pool = LR_Selection(G, J);

    %//////////////////////////////////////////////////////////////////////
    %/// CROSSOVER
    G_crossed = Crossover(G_pool, CR);

    %//////////////////////////////////////////////////////////////////////
    %/// TAGUCHI METHOD   
    G_taguchi = taguchi_method(G_crossed, CR, J);
    
    %//////////////////////////////////////////////////////////////////////
    %/// MUTATION
    G_mutated = Mutate(G_taguchi, MR);
    if(size(G_mutated,2)<P)
        for i = 1:P - size(G_mutated,2)
            G_mutated = [G_mutated G_mutated(:,1)];
        end
    end
       
    %//////////////////////////////////////////////////////////////////////
    %/// REPLACEMENT (Elitist)
    F_s = Eval_Fitness(G_mutated, J);
    func_evals = func_evals + size(G_mutated,2);
    
    G_sort = [F_s.fitness' G_mutated']';
    G_sort = sortrows(G_sort',1)';
    G = G_sort(:,1:P);
    F_s.fitness = G_sort(1,1:P);
    G(1,:) = [];

    %//////////////////////////////////////////////////////////////////////
    %/// TRACK METRICS
    avg_cost(loopcount)     = mean(F_s.fitness);
    var_track(loopcount)    = var(avg_cost);
    best_cost(:,loopcount)  = F_s.best_cost;
    worst_cost(:,loopcount) = F_s.worst_cost;
    best_chr(:,loopcount)   = F_s.best_chr;
    worst_chr(:,loopcount)  = F_s.worst_chr;
    avg_chr(:,loopcount)    = mean(G');

    metric_s.avg_cost   = avg_cost;
    metric_s.var_track  = var_track;
    metric_s.best_cost  = best_cost;
    metric_s.worst_cost = worst_cost;
    metric_s.best_chr   = best_chr;
    metric_s.worst_chr  = worst_chr;
    metric_s.avg_chr    = avg_chr;

    % Increment convergence variable if population variance stalled
%     if( (loopcount > 1) & (var(avg_cost(loopcount-1:loopcount)) < epsilon) )
    if( (loopcount > min_iter)&(best_cost(:,loopcount-1) == best_cost(:,loopcount)) )
        converge = converge + 1;
    else
        converge = 0;
    end
    
    % Maximum-Iteration Limit
    if(loopcount > max_iter)
        break
    end
    
    loopcount = loopcount + 1;

end

Darwin_s.F_s      = F_s;
Darwin_s.metric_s = metric_s;

end