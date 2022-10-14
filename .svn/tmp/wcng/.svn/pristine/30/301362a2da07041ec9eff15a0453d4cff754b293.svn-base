function [DSMoutput] = DTsimulation(NTF,STF,OSR,SIM_CONFIG)

%% Initialization
global order

WRAPS   = SIM_CONFIG.runs;    % Number of runs to average
F_C     = 1/OSR;              % DSM corner frequency
f0      = SIM_CONFIG.f0;
db_down = SIM_CONFIG.db_down; % Input signal scaling
FFT_LEN = SIM_CONFIG.FFT_LEN; % Number of FFT points

% Number of Data Points
DECIM = 2*OSR;                
N     = 2*WRAPS*DECIM*FFT_LEN;

% Input Frequency
% FFT_BIN = ceil(FFT_LEN/10);  % Low-Pass FFT_BIN
FFT_BIN = ceil(f0*FFT_LEN);
FREQ    = 8*WRAPS*FFT_BIN;

% Input Signal
input_scale = db_down/3;
A = 2^(-input_scale);
x = A*cos(2*pi*FREQ/N*(0:N-1));

%% Discrete-Time Simulation

%//////////////////////////////////////////////////////////////////////////
%/// DIFFERENCE EQUATION SOLUTION

% Loop-Variable Initialization
q_in  = zeros(1,N);
y     = zeros(1,N);
gain  = zeros(1,N);
y1    = zeros(1,N); y2    = zeros(1,N);
y1_FF = zeros(1,1); y2_FF = zeros(1,1);
y1_FB = zeros(1,1); y2_FB = zeros(1,1);

F = STF.num; 
G = NTF.num;
G = G / G(1);
H = NTF.den - G;
H = [H(2:end) 0];

for n = order + 2 : N
    
    for i = 1 : order + 1
        y1_FF(i) = F(i)*x(n-i+1);
        y2_FF(i) = H(i)*y(n-i);
    end
    
    for i = 2 : order + 1
        y1_FB(i-1) = G(i)*y1(n-i+1);
        y2_FB(i-1) = G(i)*y2(n-i+1);
    end
    
    y1(n) = sum(y1_FF)-sum(y1_FB);
    y2(n) = sum(y2_FF)-sum(y2_FB);
   
%     q_in(n) = x(n) - y2(n);
    q_in(n) = y1(n) - y2(n);
    
    y(n) = sign(q_in(n));
% 	y(n) = q_in(n);
end

DSMoutput.y    = y;
DSMoutput.q_in = q_in;
DSMoutput.gain = gain;

%% Output Decimation Filter 

%//////////////////////////////////////////////////////////////////////////
%/// LOW-PASS ANTI-ALIASING FILTER

filter_name = ['LP_OSR_' num2str(OSR) '_AA.mat'];

if isempty( which(filter_name) )
    
    FIR_FFT_LEN = 2^13;
    
    % Parks-Mclellan Optimal Equiripple FIR Design
    pb_ripple = 0.1;
    sb_ripple = 200;
    f = [F_C 2*F_C];
    a = [1 0];
    dev = [(10^(pb_ripple/20)-1)/(10^(pb_ripple/20)+1)  10^(-sb_ripple/20)];
    [FIR_ORDER,f0,a0,w]=firpmord(f,a,dev);
    b = firpm(FIR_ORDER,f0,a0,w);
    
    % FIR Frequency Response
    [FIR_jw,w] = freqz(b,1,FIR_FFT_LEN);
    FIR_mag    =  db(abs(FIR_jw));
    
    % FIR Impulse Response
    figure(100); clf
    clf
    plot(w(1:FIR_FFT_LEN/2)/pi,FIR_mag(1:FIR_FFT_LEN/2))
    grid on
    title('Decimation FIR Magnitude Response (dB)')
    xlabel('Normalized Frequency (x \pi rad/sample)')
    ylabel('Magnitude (dB)')
    
    save(filter_name,'b','FIR_mag')

else
    dec_filter = open(filter_name);
    b = dec_filter.b;
end

% Anti-Aliased Output
y_filt  = conv(b,y);
y_trunc = y_filt(1024:end-1024);


%//////////////////////////////////////////////////////////////////////////
%/// OUTPUT DECIMATION

y = downsample(y_trunc,0.5*OSR);
y = y(1:WRAPS*FFT_LEN);

%% SNR-Calculation

BAND_EDGE = FFT_LEN/4;

% Chebyshev Window
win = chebwin(FFT_LEN,200)'; bin_offset = 7;
% win = hann(FFT_LEN); bin_offset = 2;
window_P = mean(win.^2);

loop_ave = zeros(1,FFT_LEN);

for i = 1:WRAPS
    
    % FFT Wrap indexing
    START = 1 + FFT_LEN*(i-1);
    
    loop_jw  = fft(y(START:START+FFT_LEN-1).*win,FFT_LEN);
    loop_ave = loop_ave + loop_jw.*conj(loop_jw);
   
end

Y_JW  = loop_ave / WRAPS;
Y_mag = 10*log10(abs(Y_JW));
shift = max(Y_mag);
Y_mag = Y_mag - shift;

% 'In-Band' Spectrum
inband_JW = Y_JW(1:BAND_EDGE);

% Signal Bins
% signal_index = find(Y_mag==max(Y_mag(bin_offset+1:floor(length(Y_mag)/2))))
signal_index = FFT_BIN + 1;
signal_JW    = inband_JW(signal_index - bin_offset : signal_index + bin_offset);
signal_P     = sum(signal_JW)/FFT_LEN/window_P;

% Noise Bins
noise_JW  = inband_JW;
noise_JW(signal_index - bin_offset : signal_index + bin_offset) = [];
noise_P   = sum(noise_JW)/FFT_LEN/window_P;

% Dynamic Range
DR_P = BAND_EDGE*max(noise_JW)/FFT_LEN/window_P;

% Signal-to-Noise Ratio
SNR_sim = 10*log10(signal_P/noise_P);

% Spurrious Free Dynamic Range
SFDR_sim = -1*(10*log10(max(noise_JW))-shift);

% Dynamic Range
DR_sim = 10*log10(signal_P/DR_P);

DSMoutput.mag  = Y_mag;
DSMoutput.SNR  = SNR_sim;
DSMoutput.DR   = DR_sim;
DSMoutput.SFDR = SFDR_sim;
