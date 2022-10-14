function [num] = gamma2poly( gamma )

% Usage: [num] = gamma2poly( gamma )
%
% Description: Converts a paratmeric vector (chromosome) 'gamma' to a nth
% order rational function characterized by the polyomial coefficients 'num'
% and 'den'

% Vector indices
gamma_i = 1;
poly_i  = 1;

% Parameters
L = length(gamma);  % Length of Chromosome
M = mod(L-1,2);   % # of 1st-order sections
N = floor((L-1)/2);     % # of 2nd-order sections
order = 2*N + M;

% 1st-Order Transfer Function
if( M==1 )
    num(1,:) = [1 gamma(1) 0];
    gamma_i = 2;
    poly_i  = 2;
end

% 2nd-Order Transfer Fucntion Zeros
for i = poly_i:N+M    
    num(i,:) = [1 gamma(gamma_i) gamma(gamma_i+1)];
    gamma_i  = gamma_i + 2;
end

temp_num = [];

for i = 1:N+M
    if (isempty(temp_num))
        temp_num = num(i,:);
    else
        temp_num = conv(temp_num,num(i,:));
    end
end

k   = gamma(M+2*N+1);
num = k*temp_num;

end