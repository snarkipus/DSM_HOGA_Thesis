function J = HP_obj(theta)

%% nth-Order High-Pass NTF Objective Function
%  //////////////////////////////////////////////////////////////////////////
% ///
% ///       Usage: J = HP_obj(theta,OSR)
% ///
% ///   Arguments: [int] OSR:   DSM Oversampling Rate
% ///              [vec] theta: paramater vector (chromosome)
% ///
% ///     Returns: [dbl] J:     Cost Evaluation
% ///
% /////////////////////////////////////////////////////////////////////////

%% File Information
%  //////////////////////////////////////////////////////////////////////////
% ///
% /// File Name: HP_obj.m
% ///
% /// Description:
% /// M-File which evaluates the cost function for a high-pass NTF.
% ///
% /// Author: M. Jackson
% ///
% /// Changelog:
% ///
% ///    [DATE]      [VERSION]
% ///
% ///    09/04/07    0.1
% ///    generalized n'th order objective function adapted from
% ///    HP4_obj.m (static 4th order function)
% ///
% ///    09/05/07    0.11
% ///    - added penalties (removed)
% ///       - 2nd-order poles on the real axis
% ///       - pole-zero cancellation
% ///
% ///    11/12/07    0.12
% ///    - Modified cost function for passband LMS and stopband MSE
% ///
% /////////////////////////////////////////////////////////////////////////

global OSR DR

% Vector indices
theta_i = 1;
poly_i  = 1;

% Parameters
L = length(theta);  % Length of Chromosome
M = mod(L-1,4)/2;   % # of 1st-order sections
N = floor(L/4);     % # of 2nd-order sections
order = 2*N + M;

[NTFnum,NTFden] = theta2poly(theta);

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
w_NTFpass = linspace(2*pi/OSR,pi,16);

% Evaluate the NTF over [w_sb w_pb]
NTFh_sb = freqz(NTFnum,NTFden,w_NTFstop);
NTFh_pb = freqz(NTFnum,NTFden,w_NTFpass);

%% Frequency Domain Analysis

% Stop-Band Region
stop_mag  = abs(NTFh_sb); % Stop-Band Magnitude
stop_db   = db(stop_mag);

% Pass-Band Region
pass_mag  = abs(NTFh_pb); % Pass-Band Magnitude
pass_db   = db(pass_mag);

%% Objective Function

% PENALTY: NTF Dynamic Range (DR) 
if( db(max(stop_mag)) > -DR )
    J = 1e2;
    return
end

% Objective Function Weights
% alpha = 0.5;
alpha = 0;

%   [ 'Pass-Band LMS' 'Stop-Band MSE' ]
v = [      alpha         (1-alpha)    ];

% J1 = v(1)*( norm(stop_db,1)    + norm((pass_db),1)  ); % 1-Norm 
% J2 = v(2)*( norm(abs(stop_db)) + norm(abs(pass_db)) ); % 2-Norm 


J1 = v(1)*( norm(stop_mag,1)    + norm((1-pass_mag),1)  ); % 1-Norm 
J2 = v(2)*( norm(abs(stop_mag)) + norm(abs(1-pass_mag)) ); % 2-Norm 

% J1 = v(1)*norm(1-pass_mag);         % Pass-Band LMS Optimization
% J2 = v(2)*norm(abs(stop_mag).^2,1); % Stop-Band MSE Optimization



J = J1 + J2;

end
