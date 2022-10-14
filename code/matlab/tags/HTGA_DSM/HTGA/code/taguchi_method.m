function G_taguchi = taguchi_method(G_crossed, CR, J)

%//////////////////////////////////////////////////////////////////////////
%///
%/// G_taguchi= taguchi_method(G_crossed, CR)
%///
%/// Description:
%/// Intelligently create offspring via the Taguchi Method
%/// [Tsai et al,2007]
%///
%/// Arguments:     [mat]   G_crossed:  Generation Matrix (post crossover)
%///                [dbl]   CR:         Crossover Rate (constant)
%///
%/// Returns:       [mat]   G_taguchi:  New Generation Matrix
%///
%//////////////////////////////////////////////////////////////////////////

% Population Information
P = size(G_crossed,2);  % Number of Chromosomes
N = size(G_crossed,1);  % Number of Traits per Chromosome

% Generate 2-Level Orthogonal Array
OA = generateOA(N);

loopcount = 1;
G_taguchi = [];

while(loopcount < floor(0.5*P*CR))

    % Randomly select 'Mom'
    pick_Mommy = floor(P*rand())+1;
    pick_Daddy = pick_Mommy;

    % Ensure that 'Dad' is different than 'Mom'
    while( pick_Mommy == pick_Daddy)
        pick_Daddy = floor(P*rand())+1;
    end

    Mommy = G_crossed(:,pick_Mommy); % Level 1
    Daddy = G_crossed(:,pick_Daddy); % Level 2
        
    % Create Experiment Matrix
    for i=1:size(OA,2)
        for j=1:size(OA,1)

            if( OA(j,i) == 1 )
                exp_m(j,i) = Mommy(i);
            else
                exp_m(j,i) = Daddy(i);
            end

        end
    end

    % Calculate results for 'n' experiments
    for i=1:size(exp_m,1)
        n(i) = 1/J(exp_m(i,:))^2;
    end

    % Calculate the effects of the experimental factors
    for i=1:size(OA,2) 
        for j=1:size(OA,1) 

            if( exp_m(j,i) == Mommy(i) )
                E1(j,i) = n(j);
                E2(j,i) = 0;
            else
                E1(j,i) = 0;
                E2(j,i) = n(j);
            end

        end
    end

    E1 = sum(E1);
    E2 = sum(E2);

    % Generate Optimal Chromosome
    for i=1:length(E1)
        if( E1(i) > E2(i) )
            progeny(i) = Mommy(i);
        else
            progeny(i) = Daddy(i);
        end
    end
    
    % Store Chromosome in the pool
    G_taguchi = [G_taguchi progeny'];

    loopcount = loopcount + 1;
    
end

% Add the Taguchi offspring to the population
G_taguchi = [G_crossed G_taguchi];

end


