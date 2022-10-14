function J = HP_obj(theta)

%% nth-Order High-Pass NTF Objective Function
%  //////////////////////////////////////////////////////////////////////////
% ///
% ///       Usage: J = HP_obj(theta,OSR)
% ///
% ///   Arguments: [vec] theta: paramater vector (chromosome)
% ///
% ///     Returns: [dbl] J: Cost Evaluation
% ///
% /////////////////////////////////////////////////////////////////////////

global OSR DR

% Parameters
L = length(theta);  % Length of Chromosome
M = mod(L-1,4)/2;   % # of 1st-order sections
N = floor(L/4);     % # of 2nd-order sections
order = 2*N + M;

[NTFnum,NTFden] = theta2poly(theta);

% avoid null numerator issues
if sum(NTFnum) == 0
    J = 1e2;
    return
end

%% Stability Check

% All poles must be contained by the unit circle.
poles = roots(NTFden);

for i=1:order
    
    % Poles within the unit circle
    if( abs(poles(i)) > 1 )
        J = 1e3;
        return
    end

end

w_NTFstop = linspace(0,pi/OSR,32);
w_NTFpass = linspace(1.2*pi/OSR,pi,32);

% Evaluate the NTF over [w_sb w_pb]
NTFh_sb = freqz(NTFnum,NTFden,w_NTFstop);
NTFh_pb = freqz(NTFnum,NTFden,w_NTFpass);

%% Frequency Domain Analysis

% Stop-Band Region
stop_mag  = abs(NTFh_sb); % Stop-Band Magnitude

% Pass-Band Region
pass_mag  = abs(NTFh_pb); % Pass-Band Magnitude

%% Objective Function

% PENALTY: NTF Dynamic Range (DR) 
if( db(max(stop_mag)) > -DR )
    J = 1e2;
    return
end

% Objective Function Weights
alpha = 0.0;
beta  = 0.1;

J1 =          alpha*( norm(stop_mag,2)   ); % Stopband 2-Norm (SNR) 
J2 =           beta*( norm(stop_mag,inf) ); % Stopband Infinity Norm (DR) 
J3 = (1-alpha-beta)*( norm(1-pass_mag,2) ); % Passband Infinity Norm 

J = J1 + J2 + J3;

end
