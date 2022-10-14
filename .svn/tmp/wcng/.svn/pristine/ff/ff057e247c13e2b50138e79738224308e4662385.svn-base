clear all
close all

% Globals
global func_evals;

num_loops = 50;

disp(' ')
disp('Starting New Run...')
disp(' ')

tic

% GA Coefficients and Settings
config_s.P            =    200; % Population Size (# of chromosomes)
config_s.CR           =    0.1; % Cross-Over Rate
config_s.MR           =   0.02; % Mutation Rate
config_s.max_iter     = 100000; % Maximum Number of Iterations
config_s.min_iter     =   1000; % Minimum Number of Iterations
config_s.converge_min =     50; % Minimum Number of Iterations for Convergence
config_s.epsilon      =  1e-20; % Misadjustment threshold

% Objective Function
test_s.type = 'rose';
test_s.N    =     2;
tic

idx_len = realmax;

for i = 1 : num_loops
    
    Results(i) = HOGA(config_s,test_s);
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
    %     worst_cost(i,:) = track{3,i}(1:min_length);
    
end

if num_loops > 1
    plot_avg_cost  = mean(avg_cost);
    plot_avg_best  = mean(best_cost);
else 
    plot_avg_cost = avg_cost;
    plot_avg_best = best_cost;
end

% plot_avg_worst = mean(worst_cost);

disp('Completed')
disp(' ')
% disp(['Best Chromosome: ' num2str(avg_best_chr)])
toc


% 3-D Peaks() Plot
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



figure(2); clf
plot(plot_avg_cost,'LineWidth',2)
hold on
plot(plot_avg_best,'r','LineWidth',2)
%     plot(plot_avg_worst,'g','LineWidth',2)
axis square
legend('Average Chromosome','Best Chromosome')
xlabel('Generations (i)')
ylabel('Mean Cost')
% axis([0 length(plot_avg_cost) -Inf Inf])
axis square tight
% title(['P_c= ' num2str(config_s.CR),' ',...
%        'P_m= ' num2str(config_s.MR),' ',...
%        'Runs= ' num2str(num_loops)])
% hold off
% else
%     figure(2); clf
%     semilogy(plot_avg_cost,'LineWidth',2)
%     hold on
%     semilogy(plot_avg_best,'r','LineWidth',2)
% %     semilogy(plot_avg_worst,'g','LineWidth',2)
%     axis square
% %     legend('Average','Best','Worst')
%     xlabel('Generations (i)')
%     ylabel('Average Cost')
%     title({['HOGA Learning Curve: ' test_s.type ' Function'];
%            ['Runs= ' num2str(num_loops) ,...
%             ' P_c= ' num2str(config_s.CR),...
%             ' P_m= ' num2str(config_s.MR)]})
hold off
% end

run_num = 1;
good_plot = 'N';

figure(3);
while( ~strcmp(good_plot,'Y') )


if strcmp(test_s.type,'peaks')
    
    clf
    
    %Top-View (Steepest Descent)
    Z = peaks(100);
    X = linspace(-3,3,size(Z,2));
    Y = linspace(-3,3,size(Z,1));
    contour(X,Y,Z)
    hold on
    scatter(Results(run_num).metric_s.avg_chr(1,:),... % x-coordinate
        Results(run_num).metric_s.avg_chr(2,:),... % y-coordinate
        'Marker','d')
    % axis([-0.6 0.6 -2.5 -0.5])
    axis square
    xlabel('X')
    ylabel('Y')
    title(['Parametric Contour Plot for Run ' num2str(run_num)])
    
else
    
    if strcmp(test_s.type,'rast')
        fcn = @rastriginsfcn;
    else
        fcn = @dejong2fcn;
    end
    
    clf
    
    range = [-2,2;-2,2];
    
    pts = 200;
    span = diff(range')/(pts - 1);
    x = range(1,1): span(1) : range(1,2);
    y = range(2,1): span(2) : range(2,2);
    
    pop = zeros(pts * pts,2);
    k = 1;
    for i = 1:pts
        for j = 1:pts
            pop(k,:) = [x(i),y(j)];
            k = k + 1;
        end
    end
    
    values = feval(fcn,pop);
    values = reshape(values,pts,pts);
    
    contour(x,y,log(values))
    hold on
    title(['Parametric Contour Plot for Run ' num2str(run_num)])
    scatter(Results(run_num).metric_s.avg_chr(1,:),... % x-coordinate
        Results(run_num).metric_s.avg_chr(2,:),... % y-coordinate
        'MarkerEdgeColor',[0.4784 0.06275 0.8941],'Marker','^')
    grid on
    %     axis([-1 1 -1 1])
    
end

    disp('Press any key to continue...')
    pause

    if run_num < num_loops
        run_num = run_num + 1;
    else
        run_num = 1;
    end
end

clf
surf(x,y,values,'MeshStyle','row','EdgeAlpha',0.4)
% shading interp
% light
% lighting phong
axis tight
axis square
hold on
contour(x,y,values)
rotate3d
view(37,60)
% %
% reply = 'N';
% run_num = 1;
% 
% while(strcmp(reply,'Y')==0)
% 
%     figure(6); clf
%     contour(x,y,values)
%     hold on
%     title(['Parametric Contour Plot for Run ' num2str(run_num)])
%     scatter(Results(run_num).metric_s.avg_chr(1,:),... % x-coordinate
%             Results(run_num).metric_s.avg_chr(2,:),... % y-coordinate
%             'MarkerEdgeColor',[0.4784 0.06275 0.8941],'Marker','^')
% %     axis([-1 1 -1 1])
% 
% % Z = peaks(100);
% % X = linspace(-3,3,size(Z,2));
% % Y = linspace(-3,3,size(Z,1));
% % contour(X,Y,Z)
% % hold on
% % scatter(Results(run_num).metric_s.avg_chr(1,:),... % x-coordinate
% %         Results(run_num).metric_s.avg_chr(2,:),... % y-coordinate
% %         'Marker','d')
% % axis([-0.6 0.6 -2.5 -0.5])
% axis square
% xlabel('X')
% ylabel('Y')
% title('Average Chromosome Value')
% 
%     reply = input('Happy Plot? [N]:','S');
%         if isempty('S')
%                reply = 'N'
%         else
%             reply = upper('S');
%         end
% 
%         if run_num < num_loops
%             run_num = run_num + 1;
%         else
%             run_num = 1
%         end
% 
% 
% end
        