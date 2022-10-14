function Darwin_s = HOGA(config_s)

%% Hybrid Orthogonal Genetic Algorithm (HOGA) 
%  //////////////////////////////////////////////////////////////////////////
% ///
% ///       Usage: Darwin_s = HOGA(config_s, test_s)
% ///
% ///   Arguments: [struct] test_s   : Objective Function
% ///              [struct] config_s : HOGA configuration
% ///
% ///     Returns: [struct] Darwin_s : Evolved Results 
% ///
% /// Definitions:
% ///
% ///    <Darwin_s>      [struct]    F_s: Final results
% ///     |              [struct]    metric_s: Tracked results
% ///     |
% ///     |-------<F_s>  [vec] best_chr
% ///     |              [dbl] best_cost
% ///     |              [vec] worst_chr
% ///     |              [dbl] worst_cost
% ///     |              [dbl] fitness
% ///     | 
% ///     |- <metric_s>  [vec] avg_cost 
% ///                    [vec] var_track (fitness variance)
% ///                    [vec] best_cost
% ///                    [vec] worst_cost
% ///                    [vec] best_chr
% ///                    [vec] worst_chr
% ///                    [vec] avg_chr
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
% ///    02/04/09   0.51
% ///     -  Fixed comments and cleaned up code 
% ///     -  Removed elitist replication block
% ///
% //////////////////////////////////////////////////////////////////////////

%%  GA Parameter Initialization

% HOGA Configuration
P  = config_s.P;
CR = config_s.CR; 
MR = config_s.MR;

max_iter     = config_s.max_iter;
min_iter     = config_s.min_iter;
converge_min = config_s.converge_min;
epsilon      = config_s.epsilon;

%% INITIALIZATION
converge  = 0;
loopcount = 1;
best_cost = realmax;

[J,G] = config_s.init_fcn(P);
F_s   = Eval_Fitness(G,J);
    
%% MAIN LOOP
while(converge < converge_min)
    
    %% LINEAR-RANKING SELECTION
    G_pool = LR_Selection(G, F_s);

    %% CROSSOVER
    G_crossed = Crossover(G_pool, CR);

    %% HYBRID ORTHOGONAL CROSSOVER VIA TAGUCHI METHOD   
    G_taguchi = taguchi_method(G_crossed, CR, J);
    
    %% MUTATION
    G_mutated = Mutate(G_taguchi, MR);
    
    %% EVALUATION
    F_s = Eval_Fitness(G_mutated, J);   

    %% DATA TRACKING
    metric_s.avg_cost(loopcount)     = mean(F_s.fitness);
    metric_s.var_track(loopcount)    = var(metric_s.avg_cost(loopcount));
    metric_s.best_cost(:,loopcount)  = F_s.best_cost;
    metric_s.worst_cost(:,loopcount) = F_s.worst_cost;
    metric_s.best_chr(:,loopcount)   = F_s.best_chr;
    metric_s.worst_chr(:,loopcount)  = F_s.worst_chr;
    metric_s.avg_chr(:,loopcount)    = mean(G_mutated');

    %% CONVERGENCE
    if (loopcount > min_iter)
        if(F_s.best_cost == min(metric_s.best_cost))
            converge = converge + 1;
        else
            converge = 0;
        end
    else
        converge = 0;
    end   
    
    % Maximum-Iteration Limit
    if(loopcount > max_iter); break; end
    
    loopcount = loopcount + 1;
    G = G_mutated;

end

Darwin_s.F_s      = F_s;
Darwin_s.metric_s = metric_s;

end