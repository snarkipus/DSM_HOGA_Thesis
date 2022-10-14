function [J, G_init] =  Setup_Obj_Fcn(P, N, type)

%% Setup_Obj_Fcn
%  //////////////////////////////////////////////////////////////////////////
% ///
% ///       Usage: [J, G_init] = Setup_Obj_Fcn(P, N, type)
% ///
% ///   Arguments: [dbl]          P: Population Size
% ///              [dbl]          N: Chromosome Length (Order)
% ///              [str]       type: objective function type string
% ///
% ///     Returns: [fcn_hdl]      J: Objective Function Handle
% ///              [mat]     G_init: Initialized Population
% ///
% //////////////////////////////////////////////////////////////////////////

% Variable Preallocation and Initialization
G_init     = zeros(N,P); 
theta      = zeros(1,N);
theta_min  = zeros(1,N); 
theta_max  = zeros(1,N);

switch upper(type)

    case 'PEAKS'
        %//////////////////////////////////////////////////////////////////////////
        %/// Peaks(x,y) - Optimal Solution: -6.5466@[0.25,-1.625]'ish
        J = @(vec_in) peaks(vec_in(1),vec_in(2));

        %Constraints: -3 <= (x,y) <= 3
        theta_min(1) = -3; theta_max(1) = 3;
        theta_min(2) = -3; theta_max(2) = 3;

    case 'RAST'
        %//////////////////////////////////////////////////////////////////////////
        %/// Rastigrin's Function - Optimal Solution: 0@[0,0]

        J = @(vec_in) N*10 + sum(vec_in.^2 - 10*cos(2*pi*vec_in));
        
        for i = 1:N
            theta_min(i) = -5.12; theta_max(i) = 5.12;
        end
        
    case 'ROSE'
        %//////////////////////////////////////////////////////////////////////////
        %/// Rosenbrock Function - Optimal Solution: [1,1]       
        J = @(vec_in) 100*(vec_in(1)^2-vec_in(2))^2+(1-vec_in(1))^2;
        
        %Constraints: -3 <= (x,y) <= 3
        theta_min(1) = -3; theta_max(1) = 3;
        theta_min(2) = -3; theta_max(2) = 3;
        
        
    case 'PARABOLA'
        %//////////////////////////////////////////////////////////////////////////
        %/// f_11 : Jmin = 0
        J = @(vec_in)   sum(vec_in.^2);
        
        for i = 1:N
            theta_min(i) = -100; theta_max(i) = 100;
        end
end

% % %//////////////////////////////////////////////////////////////////////////
% %/// f_1
% J = @(vec_in) sum(-vec_in.*sin(sqrt(abs(vec_in))));
%
% % Constraints: [-500,500]^N
% for i = 1:N
%     theta_min(i) = -500; theta_max(i) = 500;
% end

%//////////////////////////////////////////////////////////////////////////
%/// f_9
% J = @(vec_in)   (1/N)*sum(vec_in.^4-16*vec_in.^2+5*vec_in);
%
% % Constraints: [-5,5]^N
% for i = 1:N
%     theta_min(i) = -5; theta_max(i) = 5;
% end

% %//////////////////////////////////////////////////////////////////////////
% %/// f_11 : Jmin = 0
% J = @(vec_in)   sum(vec_in.^2);

% Generate the Initial Random Population
for i = 1:P
    for j = 1:N
        theta(j) = theta_min(j) + rand()*(theta_max(j) - theta_min(j));
    end
    G_init(:,i) = theta';
end

end