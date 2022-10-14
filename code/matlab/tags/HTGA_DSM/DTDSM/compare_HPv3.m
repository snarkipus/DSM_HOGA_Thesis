clear
close all

%% Initialization


timestamp = clock;
disp(' ')
disp(['START TIME >> ' num2str(timestamp(4)) ':' num2str(timestamp(5))]) 
tic

order   =    5;
OSR     =   128;
runs    =    2;
FFT_LEN = 2^10;

%% DSM System Design

% Delsig Toolbox
[dsNTF,dsSTF,DR] = delsigLP(order,OSR);

% Chebyshev Design
[chNTF,chSTF] = chebyLP(order,OSR,DR+3);

%% DSM System Frequency Reponse

[dsNTF_jw,w] = freqz(dsNTF.num,dsNTF.den,FFT_LEN);
dsSTF_jw     = freqz(dsSTF.num,dsSTF.den,FFT_LEN);
dsNTF_mag    = 20*log10(abs(dsNTF_jw));
dsSTF_mag    = 20*log10(abs(dsSTF_jw));

chNTF_jw     = freqz(chNTF.num,chNTF.den,FFT_LEN);
chSTF_jw     = freqz(chSTF.num,chSTF.den,FFT_LEN);
chNTF_mag    = 20*log10(abs(chNTF_jw));
chSTF_mag    = 20*log10(abs(chSTF_jw));

%% DSM Magnitude Response Plots

% NTF Magnitude Repsonses
figure(1);clf
plot(w./pi,dsNTF_mag,'--b',w./pi,chNTF_mag,'b')
title({'\Delta\SigmaM NTF Magnitude Response Comparison';...
      ['Order:' num2str(order)...
       '  OSR:' num2str(OSR)]})
xlabel('Normalized Frequency (x \pi rad/sample)')
ylabel('Magnitude (dB)')
legend('DelSig NTF','ChebyII NTF','Location','southeast')
axis([0 .2 min(min(dsNTF_mag),min(chNTF_mag)) 10])
grid on

% STF Magnitude Responses
figure(2);clf
plot(w./pi,dsSTF_mag,'--r',w./pi,chSTF_mag,'r')
title({'\Delta\SigmaM STF Magnitude Response Comparison';...
      ['Order:' num2str(order)...
       '  OSR:' num2str(OSR)]})
xlabel('Normalized Frequency (x \pi rad/sample)')
ylabel('Magnitude (dB)')
legend('DelSig STF','ChebyII STF','Location','southeast')
axis([0 .5 min(chSTF_mag) 10])
grid on

%% DSM Simulation [DTsimulation.m]

% FFT_LEN = 2^13;

dsOutput = DTsimulation(dsNTF,dsSTF,OSR,runs,FFT_LEN);
chOutput = DTsimulation(chNTF,chSTF,OSR,runs,FFT_LEN);

table1 = ...
{'SIMULATION', 'RESULTS' , ['ORDER: ' num2str(order)], ['OSR: ' num2str(OSR)];...
 '**********' ,'*******' , '********', '********';... 
 'TYPE'   , 'SNR[dB]'   ,'DR[dB]'    , 'FLOOR[dB]';...
 'DelSig' , dsOutput.SNR, dsOutput.DR, dsOutput.SFDR;...
 'ChebyII', chOutput.SNR, chOutput.DR, chOutput.SFDR};

disp(' ')
disp(table1)

%% Undecimated Output Spectrum

Y_jw  = fft(chOutput.y(1024:1024+FFT_LEN-1).*hann(FFT_LEN));
Y_mag = db(abs(Y_jw));
w     = 0:FFT_LEN^-1:2*pi;

figure(3); clf
plot(w(1:FFT_LEN/2)/pi,Y_mag(1:FFT_LEN/2))
grid on
title({'Undecimated \Delta\SigmaM Output Spectrum';...
      ['Order:' num2str(order)...
       ' OSR:' num2str(OSR)]})
xlabel('Normalized Frequency (x \pi rad/sample)')
ylabel('Magnitude (dB)')

%% Output Spectrum Plot

figure(4); clf
n = (1:length(dsOutput.mag)).*(4/(FFT_LEN*OSR));
plot(n,dsOutput.mag,n,chOutput.mag)
title({'\Delta\SigmaM Output Spectrum Comparison';...
      ['Order:' num2str(order)...
       ' OSR:' num2str(OSR)...
       ' Runs:' num2str(runs)...
       ' FFT:' num2str(FFT_LEN)]})
axis([0 2/OSR min(dsOutput.mag) 10])
legend('DelSig','ChebyII')
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency (x\pi rad/sample)')

toc

%% Variable Storage

header = [ num2str(order) '_' num2str(OSR) '_LP' ];
save header;


