% ICSENG HTGA SCRIPT

clear all
close all

global func_evals;
func_evals = 0;

num_loops = 2;

disp(' ')
disp('Starting New Run...')
disp(' ')

tic

% GA Coefficients
config_s.P            =   200; % Population Size (# of chromosomes)
config_s.CR           =   0.8; % Cross-Over Rate
config_s.MR           =  0.5; % Mutation Rate
config_s.max_iter     = 100000; % Maximum Number of Iterations
config_s.min_iter     = 10000; % Minimum Number of Iterations
config_s.converge_min =    50; % Minimum Number of Iterations for Convergence
config_s.epsilon      = 1e-20; % Misadjustment threshold

% Objective Function
test_s.type = 'rose';
test_s.N    =     30;
tic

idx_len = realmax;

for i = 1 : num_loops
    
    Results(i)    = HOGA(config_s,test_s);

    track{1,i} = Results(i).metric_s.avg_cost;
    track{2,i} = Results(i).metric_s.best_cost;
    track{3,i} = Results(i).metric_s.worst_cost;
    track{4,i} = length(track{1,i});
    track{5,i} = Results(i).metric_s.best_chr;
    
end

min_length = min([track{4,:}]);

% [best_val,best_index] = min(track{2,:}');

% best_chr = track{5,best_index};

for i = 1 : num_loops
    
    avg_cost(i,:)   = track{1,i}(1:min_length);
    best_cost(i,:)  = track{2,i}(1:min_length);
    worst_cost(i,:) = track{3,i}(1:min_length);

end

plot_avg_cost  = mean(avg_cost);
plot_avg_best  = mean(best_cost);
plot_avg_worst = mean(worst_cost);

disp('Completed')
disp(' ')
% disp(['Best Chromosome: ' num2str(avg_best_chr)])
toc


% % 3-D Peaks() Plot
% figure(1); clf
% [X,Y,Z] = peaks(200);
% surfc(X,Y,Z,'MeshStyle','row','EdgeAlpha',0.4);
% % axis([-3 3 -3 3 -10 8])
% axis tight
% axis square
% view([-50 15])
% box on
% grid off
% title('Peaks(x,y)')

if strcmp(test_s.type,'peaks')
    
    figure(2); clf
    plot(plot_avg_cost,'LineWidth',2)
    hold on
    plot(plot_avg_best,'r','LineWidth',2)
    plot(plot_avg_worst,'g','LineWidth',2)
    axis square
    legend('Average','Best','Worst')
    xlabel('Generations (i)')
    ylabel('Average Cost')
    title({['HOGA Learning Curve: ' test_s.type ' Function'];
           ['Runs= ' num2str(num_loops) ,...
            ' P_c= ' num2str(config_s.CR),...
            ' P_m= ' num2str(config_s.MR)]})
    hold off
else
    figure(2); clf
    semilogy(plot_avg_cost,'LineWidth',2)
    hold on
    semilogy(plot_avg_best,'r','LineWidth',2)
    semilogy(plot_avg_worst,'g','LineWidth',2)
    axis square
    legend('Average','Best','Worst')
    xlabel('Generations (i)')
    ylabel('Average Cost')
    title({['HOGA Learning Curve: ' test_s.type ' Function'];
           ['Runs= ' num2str(num_loops) ,...
            ' P_c= ' num2str(config_s.CR),...
            ' P_m= ' num2str(config_s.MR)]})
    hold off
end

% random_run = floor(rand()*num_loops);
% 
% % Top-View (Steepest Descent)
% subplot(122)
% Z = peaks(100);
% X = linspace(-3,3,size(Z,2));
% Y = linspace(-3,3,size(Z,1));
% contour(X,Y,Z)
% hold on
% scatter(Results(random_run).metric_s.avg_chr(1,:),... % x-coordinate
%         Results(random_run).metric_s.avg_chr(2,:),... % y-coordinate
%         'Marker','d')
% axis([-0.6 0.6 -2.5 -0.5])
% axis square
% xlabel('X')
% ylabel('Y')
% title('Average Chromosome Value')
