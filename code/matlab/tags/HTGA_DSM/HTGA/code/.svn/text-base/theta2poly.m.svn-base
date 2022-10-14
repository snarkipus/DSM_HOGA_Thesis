function [num,den] = theta2poly( theta )

% Usage: [num,den] = theta2poly( theta )
%
% Description: Converts a paratmeric vector (chromosome) 'theta' to a nth
% order rational function characterized by the polyomial coefficients 'num'
% and 'den'

% Vector indices
theta_i = 1;
poly_i  = 1;

% Parameters
L = length(theta);  % Length of Chromosome
M = mod(L-1,4)/2;   % # of 1st-order sections
N = floor(L/4);     % # of 2nd-order sections

% 1st-Order Transfer Function
if( M==1 )
    num(1,:) = [1 theta(1) 0];
    den(1,:) = [1 theta(2) 0];
    theta_i = 3;
    poly_i  = 2;
end

% 2nd-Order Transfer Fucntion Zeros
for i = poly_i:N+M    
    num(i,:) = [1 theta(theta_i) theta(theta_i+1)];
    theta_i  = theta_i + 2;
end

% 2nd-Order Transfer Function Poles
for i = poly_i:N+M    
    den(i,:) = [1 theta(theta_i) theta(theta_i+1)];
    theta_i = theta_i + 2;    
end

temp_num = [];
temp_den = [];

for i = 1:N+M
    if (isempty(temp_num))
        temp_num = num(i,:);
        temp_den = den(i,:);
    else
        temp_num = conv(temp_num,num(i,:));
        temp_den = conv(temp_den,den(i,:));
    end
end

if(M==1)
    temp_num(end) = [];
    temp_den(end) = [];
end

k   = theta(2*(M+2*N)+1);
num = k*temp_num;
den = temp_den;

end