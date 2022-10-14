function J = LPSTF_obj(gamma)

%% nth-Order Low-Pass STF Objective Function
%  //////////////////////////////////////////////////////////////////////////
% ///
% ///       Usage: J = LPSTF_obj(gamma,OSR)
% ///
% ///   Arguments: [vec] gamma: paramater vector (chromosome)
% ///
% ///     Returns: [dbl] J:     Cost Evaluation
% ///
% /////////////////////////////////////////////////////////////////////////

%% File Information
%  //////////////////////////////////////////////////////////////////////////
% ///
% /// File Name: LPSTF_obj.m
% ///
% /// Description:
% /// M-File which evaluates the cost function for a low-pass STF.
% ///
% /// Author: M. Jackson
% ///
% /// Changelog:
% ///
% ///    [DATE]      [VERSION]
% ///
% ///    09/09/07    0.1 (initial work)
% ///
% /////////////////////////////////////////////////////////////////////////

global OSR DR myNTFden

% Vector indices
theta_i = 1;
poly_i  = 1;

% Parameters
L = length(gamma);  % Length of Chromosome
M = mod(L-1,2)/2;   % # of 1st-order sections
N = floor(L/2);     % # of 2nd-order sections
order = 2*N + M;

[STFnum] = gamma2poly(gamma);

w_STFpass = linspace(0,pi/OSR,32);
w_STFstop = linspace(2*pi/OSR,pi,16);

% Evaluate the NTF over [w_sb w_pb]
STFh_pb = freqz(STFnum,myNTFden,w_STFpass);
STFh_sb = freqz(STFnum,myNTFden,w_STFstop);

%% Frequency Domain Analysis

% Pass-Band Region
pass_mag  =   abs(STFh_pb); % Pass-Band Magnitude

% Stop-Band Region
stop_mag  =   abs(STFh_sb); % Stop-Band Magnitude

%% Objective Function

% Objective Function Weights

alpha = 1;

%   [ '1-norm' '2-norm' ]
v = [  alpha  (1-alpha) ];

% PENALTY: NTF Dynamic Range (DR) 
if( db(stop_mag(end)) > -0.5*DR )
    J = 1e2;
    return
end

% PENALTY: STF Peaking
if( max( max(db(pass_mag)),max(db(stop_mag)) ) > 1.1)
    J = 1e3;
    return
end

J1 = v(1)*(norm(stop_mag,1)    + norm(1-pass_mag,1));    % 1-Norm 
J2 = v(2)*(norm(abs(stop_mag)) + norm(abs(1-pass_mag))); % 2-Norm 

% J1 = v(1)*(norm(1-pass_mag,1));  % 1-Norm 
% J2 = v(2)*(norm(abs(stop_mag))); % 2-Norm 

J = J1 + J2;

end
