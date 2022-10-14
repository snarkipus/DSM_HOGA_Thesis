function J = LPSTF_obj(gamma)

%% nth-Order Low-Pass STF Objective Function
%  //////////////////////////////////////////////////////////////////////////
% ///
% ///       Usage: J = LPSTF_obj(gamma,OSR)
% ///
% ///   Arguments: [vec] gamma: paramater vector (chromosome)
% ///
% ///     Returns: [dbl] J: Cost Evaluation
% ///
% /////////////////////////////////////////////////////////////////////////

global OSR DR myNTFden

[STFnum] = gamma2poly(gamma);

% avoid null numerator issues
if sum(STFnum) == 0
    J = 1e2;
    return
end

w_STFpass = linspace(0,pi/OSR,32);
w_STFstop = linspace(2*pi/OSR,pi/2,32);

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

% PENALTY: NTF Dynamic Range (DR) 
if( db(stop_mag(end-1)) > -0.9*DR )
    J = 1e3;
    return
end

% PENALTY: STF Peaking
if( max( max(db(pass_mag)),max(db(stop_mag)) ) > 1.1)
    J = 1e3;
    return
end

J1 =     alpha*( norm(1-pass_mag,1) ); % Passband 1-Norm (Maximally Flat) 
J2 = (1-alpha)*( norm(stop_mag,2)   ); % Stopband 2-Norm (SNR)

J = J1 + J2;

end
