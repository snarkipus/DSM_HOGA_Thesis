clear all
close all

set(0,'defaulttextinterpreter','none')

%% Initialization

FFT_LEN = 2^13;          % FFT Length
B       = 6;             % Number of Quantization Bits
WRAPS   = 10;         % Number of runs to average
N       = WRAPS*FFT_LEN; % Signal Length

% Input Frequency
FFT_BIN = 395;
FREQ    = WRAPS*FFT_BIN;

% Dithered Input Signal (triangular)
s      = (1-2^-B)*cos(2*pi*FREQ/N*[0:N-1]);
% dither = 2^(-B+1)*(rand(1,N)-0.5 + rand(1,N)-0.5);
dither = zeros(1,N);
x      = s + dither;

%% Quantization

% A`la Stubberud
[xq,e] = quant(x + 2^-B,B);
xq     = xq - 2^-B;

y = xq;

%% SNR-Calculation

BAND_EDGE = FFT_LEN/2;

% Chebyshev Window
win = chebwin(FFT_LEN,125)'; bin_offset = 10;
% win = hann(FFT_LEN); bin_offset = 1;
window_P = mean(win.^2);

loop_ave    = zeros(1,FFT_LEN);

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
inband_JW =   Y_JW(1:BAND_EDGE);

% Signal Bins
signal_index = find(Y_mag==max(Y_mag));
signal_JW    = inband_JW(signal_index - bin_offset : signal_index + bin_offset);
% signal_JW = inband_JW(signal_index(1));
signal_P     = sum(signal_JW)/FFT_LEN/window_P

% Noise Bins
noise_JW  = inband_JW;
noise_JW(signal_index - bin_offset : signal_index + bin_offset) = [];
noise_P   = sum(noise_JW)/FFT_LEN/window_P

% Dynamic Range
DR_P = BAND_EDGE*max(noise_JW)/FFT_LEN/window_P

% Signal-to-Noise Ratio
SNR_sim = 10*log10(signal_P/noise_P);
SNRs    = sprintf('SNR   = %4.2fdB',SNR_sim);

% Noise Floor (average)
FLOOR_sim = mean(10*log10(noise_JW)-shift);
FLOORs    = sprintf('FLOOR = %4.2fdB',FLOOR_sim);

% Dynamic Range
DR_sim = 10*log10(signal_P/DR_P);
DRs    = sprintf('DR     = %4.2fdB',DR_sim);

%% Output

% Annotations
peak_index = find(Y_mag==max(Y_mag))
first_peak = peak_index(1);

front_Y = Y_mag(1:length(Y_mag)/2);
front_noise = [front_Y(1:first_peak-5) front_Y(first_peak+5:end)];
max_floor = max(front_noise);
noise_index = find(front_noise==max_floor);


figure(1);clf
plot(Y_mag)
title({'Nyquist Converter Output Spectrum';
       '(undithered sinusoidal input)'})
% title({'Nyuist Converter Ouptput Spectrum';...
%       [' Bits:' num2str(B)...
%        ' Runs:' num2str(WRAPS)...
%        ' FFT:'  num2str(FFT_LEN)]})
% text(300, -5,SNRs,  'Fontsize',12)
% text(300,-10,DRs,   'Fontsize',12)
% text(300,-15,FLOORs,'Fontsize',12)
axis([0 FFT_LEN/2 FLOOR_sim 5])
ylabel('Magnitude (dB)')
% grid on
line([peak_index(1)+100 noise_index+200 ],[max_floor max_floor],'LineWidth',1,'Color','k');
line([peak_index(1)+100 noise_index+200 ],[0 0],'LineWidth',1,'Color','k');
set(gca,'XTick',0:FFT_LEN/4:FFT_LEN/2)
set(gca,'XTickLabel',{'0','$\pi/4$','$\pi/2$'})

SQNRtheory = 6.02*B + 1.76;
SQNRdither = SQNRtheory - 4.8;

text = {'SNR',   'DR',   'FLOOR',   'Theor. SNR','Theor. SNR (dither)';...
         SNR_sim, DR_sim, FLOOR_sim, SQNRtheory,  SQNRdither};
disp(text)
