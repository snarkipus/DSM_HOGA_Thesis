function [snr,amp] = simulateSNR(ntf,OSR,amp,f0,nlev,f,k)
%[snr,amp] = simulateSNR(ntf,OSR,amp,f0=0,nlev=2,f=1/(4*OSR),k=13)
%Determine the SNR for a delta-sigma modulator by using simulations.
%The modulator is described by a noise transfer function (ntf)
%and the number of quantizer levels (nlev).
%Alternatively, ntf may be an ABCD matrix or 
%ntf can be the name of a function which takes a single argument
%representing the input signal.
%The band of interest is defined by the oversampling ratio (OSR)
%and the center frequency (f0).
%The input signal is characterized by the amp vector and the f variable.
%amp defaults to [-120 -110...-20 -15 -10 -9 -8 ... 0]dB, where 0 dB means
%a full-scale (peak value = nlev-1) sine wave.
%f is the input frequency, normalized such that 1 -> fs;
%f is rounded to an FFT bin.
%
%Using sine waves located in FFT bins, the SNR is calculated as the ratio
%of the sine wave power to the power in all in-band bins other than those
%associated with the input tone. Due to spectral smearing, the input tone
%is not allowed to lie in bins 0 or 1. The length of the FFT is 2^k.
%
% Future versions may accommodate STFs.

% Handle the input arguments
if nargin<1
    error('Insufficient arguments');
end
parameters = {'ntf';'OSR';'amp';'f0';'nlev';'f';'k'};
defaults = {NaN 64 NaN 0 2 NaN 13};
for arg_i=1:length(defaults)
    parameter = char(parameters(arg_i));
    if arg_i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
	 eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
	    eval([parameter '=defaults{arg_i};'])
    end
end
if isnan(amp)
    amp = [-120:10:-20 -15 -10:0];
end
if f0==0
    osr_mult=2;
else
    osr_mult=4;
end
if isnan(f)
    f = f0 + 1/(2*OSR*osr_mult);
end

if abs(f-f0) > 1/(OSR*osr_mult)
    fprintf(1,'Warning: the input tone is out-of-band.\n');
end

N = 2^k;
if N < 8*2*OSR	% Require at least 8 bins to be "in-band"
    fprintf(1,'Warning: Increasing k to accommodate a large oversampling ratio.\n');
    k = ceil(log2(8*2*OSR))
    N = 2^k;
end
F = round(f*N);
if F<=1
    fprintf(1,'Warning: Increasing k to accommodate a low input frequency.\n');
    % Want f*N >= 1
    k = ceil(log2(1/f))
    N = 2^k;
    F = 2;
end

Ntransient = 100;
tone = (nlev-1) * sin(2*pi*F/N*[0:(N+Ntransient-1)]);
window = .5*(1 - cos(2*pi*(0:N-1)/N) );		%Hann window
if f0==0
    % Exclude DC and its adjacent bin
    inBandBins = [3:round(N/(2*OSR))];
    F = F-2;
else
    f1 = max(round(N*(f0-1/(4*OSR))),1);
    inBandBins = [f1:round(N*(f0+1/(4*OSR)))];
    F = F-f1+1;
end

snr = zeros(size(amp));
i=1;
for A = 10.^(amp/20);
	if ~ischar(ntf)
	    v = simulateDSM(A*tone, ntf, nlev);
	else
	    v = feval(ntf, A*tone);
	end
    hwfft = fft(window.*v(1+Ntransient:N+Ntransient));
    snr(i) = calculateSNR(hwfft(inBandBins),F);
    i=i+1;
end
