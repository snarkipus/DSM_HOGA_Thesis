% Delta-Sigma Toolbox
% Version 7.1 2007.03.16
% R. Schreier, Analog Devices Inc.
%
% Consult the DSToolbox.pdf file for more complete documentation.
%
%
% Modulator synthesis and simulation
%   synthesizeNTF - Noise transfer function (NTF) synthesis.
%   clans         - Closed-loop analysis of noise shapers 
%		    (NTF synthesis for multi-bit modulators).
%   simulateDSM   - Simulate a delta-sigma modulator using either 
%		    its NTF or its state-space description.
%   simulateSNR   - Use simulateDSM to simulate a DSM with sine wave inputs 
%		    of varying amplitudes and then determine the SNR for each.
%   predictSNR    - SNR predicttion for binary modulators 
%		    (Ardalan & Paulos method)
%
% Modulator realization
%   realizeNTF	  - Compute coefficients for a particular modulator topology.
%   stuffABCD	  - Create the state-space description of a modulator given
%                   the coefficients for a particular topology.
%   scaleABCD     - Perform dynamic-range scaling.
%   mapABCD       - Convert a state-space description back to coefficients.
%
% Other functions related to delta-sigma 
%   simulateESL   - Simulate the element selection logic in a mismatch-shaping
%                   DAC.
%   designHBF     - Design multiplierless half-band filters which use the 
%                 - the Saramaki recursive filter structure.
%   findPIS       - Compute a positively-invariant set for a DSM. (The 
%                   PosInvSet sub-directory will need to be added to your PATH)
%   infnorm       - Calculate the infinity norm (maximum gain) of a TF.
%
% Demonstrations and Examples
%   dsdemo1       - Synthesize a 5th-order lowpass and an 8th-order bandpass NTF.
%   dsdemo2       - Time-domain simulation and SNR calculation.
%   dsdemo3       - Modulator realization and dynamic range scaling.
%   dsdemo4       - Audio demonstration of MOD1 and MOD2
%   dsdemo5	  - Simulate the element selection logic of a mismatch-shaping DAC.
%   dsdemo6	  - Design a hardware-efficient halfband filter.
%   dsdemo7	  - Find positively-invariant sets for second-order (and third-order) modulators.
%   dsdemo8       - Continuous-time bandpass modulator design using LC tanks.
%   dsexample1	  - Discrete-time lowpass modulator.
%   dsexample2	  - Discrete-time bandpass modulator.
%   dsexample3	  - Continuous-time bandpass modulator, including ADICE .ckt and .use files. 
%		    (ADI internal use only.)

%       Copyright (c) 1993-2007 R. Schreier
