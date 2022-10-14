function G_crossed = Crossover(G_pool, CR)

%//////////////////////////////////////////////////////////////////////////
%///
%/// G_crossed = Crossover(G_pool, CR)
%///
%/// Description:
%/// Replenish the population (post-selection) using a one-cut point
%/// operator for genetic crossover between two randomly selected progenitors.
%///
%/// Arguments:     [mat]   G_pool:     Interim Generation Matrix
%///                [dbl]   CR:         Crossover Rate (constant)
%///
%/// Returns:       [mat]   G_crossed:  New Generation Matrix
%///
%//////////////////////////////////////////////////////////////////////////

% Population Information
P = size(G_pool,2);  % Number of Chromosomes
N = size(G_pool,1);  % Number of Traits per Chromosome

Children = [];

i = 1;

pool_size = size(G_pool,2);

for j = 1 : P

    proposition = rand();

    if( proposition < CR)

        % Randomly select 'Mom'
        pick_Mommy = floor(P*rand())+1;
        pick_Daddy = pick_Mommy;

        % Ensure that 'Dad' is different than 'Mom'
        while( pick_Mommy == pick_Daddy)
            pick_Daddy = floor(P*rand())+1;
        end

        Mommy = G_pool(:,pick_Mommy);
        Daddy = G_pool(:,pick_Daddy);

        cut = floor(N*rand())+1;  % random crossover point

        % Perform Crossover
        Baby1 = [Mommy(1:cut)' Daddy(cut+1:end)']';
        Baby2 = [Daddy(1:cut)' Mommy(cut+1:end)']';

        % Have Children
        Children(:,i)   = Baby1;
        Children(:,i+1) = Baby2;

        i = i + 2;

    end

end

% Add Offspring to Population
G_crossed = [G_pool Children];

end
