function G_mutated = Mutate(G_crossed, MR)

%//////////////////////////////////////////////////////////////////////////
%///
%/// G_mutated =  Mutate(G_crossed, MR)
%///
%/// Description:
%/// Introduce genetic diversity into the population by randomly mutating
%/// chromosomes using a point-swap method.
%///
%/// Arguments:     [mat]   G_crossed:  New Generation Matrix
%///                [dbl]   MR:         Mutation Rate
%///
%/// Returns:       [mat]   G_mutated:  Mutated Generation Matrix
%///
%//////////////////////////////////////////////////////////////////////////

% Population Information
P = size(G_crossed,2);  % Number of Chromosomes
N = size(G_crossed,1);  % Number of Traits per Chromosome

G_mutated = G_crossed;

for i = 1:P

    mutagen = rand();

    if( mutagen < MR )

        idxMutate = floor(P*rand())+1;
        target    = G_mutated(:,idxMutate);

        gene_1 = floor(N*rand())+1;
        gene_2 = gene_1;

        while(gene_1 == gene_2)
            gene_2 = floor(N*rand())+1;
        end

        beta = rand();

        target_copy = target;

        target(gene_1) = (1-beta)*target_copy(gene_1) +     beta*target_copy(gene_2);
        target(gene_2) =     beta*target_copy(gene_1) + (1-beta)*target_copy(gene_2);

        G_mutated(:,idxMutate) = target;

    end

end

end