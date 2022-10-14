function [F_s, G_init] =  init_algo(P, N)

% %//////////////////////////////////////////////////////////////////////////
% %/// Peaks(x,y) - Optimal Solution: [-1.67,0.27]'ish
% J = @(vec_in) peaks(vec_in(1),vec_in(2));
% 
% %Constraints: -3 <= (x,y) <= 3
% theta_min(1) = -3; theta_max(1) = 3;
% theta_min(2) = -3; theta_max(2) = 3;

% %//////////////////////////////////////////////////////////////////////////
% %/// Rastigrin's Function - Optimal Solution: [0,0]
% J = @(vec_in)  20 + vec_in(1)^2 ...
%                   + vec_in(2)^2 ...
%                   - 10*(cos(2*pi*vec_in(1)) + cos(2*pi*vec_in(2)));
% 
% %Constraints: -5 <= (x,y) <= 5
% theta_min(1) = -5; theta_max(1) = 5;
% theta_min(2) = -5; theta_max(2) = 5;

% % %//////////////////////////////////////////////////////////////////////////
% %/// f_1
% J = @(vec_in) sum(-vec_in.*sin(sqrt(abs(vec_in))));
% 
% % Constraints: [-500,500]^N
% for i = 1:N
%     theta_min(i) = -500; theta_max(i) = 500;
% end

% %//////////////////////////////////////////////////////////////////////////
% %/// f_2
% J = @(vec_in)   sum(vec_in.^2-10*cos(2*pi*vec_in)+10);
% 
% % Constraints: [-5.12,5.12]^N
% for i = 1:N
%     theta_min(i) = -5.12^N; theta_max(i) = 5.12^N;
% end

% %//////////////////////////////////////////////////////////////////////////
% %/// f_7
% J = @(vec_in)   -sum(sin(vec_in)*sin();
%
% % Constraints: [0,PI]^N
% for i = 1:N
%     theta_min(i) = 0^N; theta_max(i) = pi^N;
% end

%//////////////////////////////////////////////////////////////////////////
%/// f_9
J = @(vec_in)   (1/N)*sum(vec_in.^4-16*vec_in.^2+5*vec_in);

% Constraints: [-5,5]^N
for i = 1:N
    theta_min(i) = -5; theta_max(i) = 5;
end

% %//////////////////////////////////////////////////////////////////////////
% %/// f_11 : Jmin = 0
% J = @(vec_in)   sum(vec_in.^2);

% Rastrigin's Function
% - Multimodal
% - High Dimension

% J = @(vec_in) N*10 + sum(vec_in.^2 - 10*cos(2*pi*vec_in));
% 
% for i = 1:N
%     theta_min(i) = -5.12; theta_max(i) = 5.12;
% end

% Assign Function Handle and Function Bounds
J_s.handle    = J;
J_s.min_bound = theta_min;
J_s.max_bound = theta_max;

% Variable Preallocation and Initialization
G_init      = zeros(N,P); % Population Matrix
F_s.fitness = zeros(1,P); % Fitness Vector
theta       = zeros(1,N); % Trait Vector
best_chr    = zeros(1,N); % Best Chromosome
best_cost   =    realmax; % Best Fitness Value
worst_chr   = zeros(1,N); % Worst Chromosome
worst_cost  =   -realmax; % Worst Fitness Value

% Process Arguments
theta_min = J_s.min_bound;
theta_max = J_s.max_bound;
J = J_s.handle;

% Generate the Initial Random Population
for i = 1:P
    for j = 1:N
        theta(j) = theta_min(j) + randn()*(theta_max(j) - theta_min(j));
    end
    G_init(:,i) = theta';
end

% Evaluate the Initial Population Fitness
F_s = Eval_Fitness(G_init, J);

% x = theta_min:0.05:theta_max;
% y = x;
% 
% for i = 1:length(x)
%     for j = 1:length(y)
%         z(i,j) = J([x(i) y(j)]);
%     end
% end
% 
% min(min(z))
% 
% figure(10)
% clf
% surfc(x,y,z,'LineStyle','none')

end