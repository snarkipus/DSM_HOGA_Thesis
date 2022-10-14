clear
close all

% GLOBAL VARIABLES
global OSR DR order myNTFden

%% Design Parameters
order =    8;               % Filter Order
OSR   =  128;               % OverSampling Rate
N     = 2^14;               % Number of Evaluation Points
fB    = ceil(N/(2*OSR));    % Noise Bandwidth (index)

makeSTF = 1;                % STF Flag

%% DelSig Toolbox NTF
ds_NTF    = synthesizeNTF(order,OSR,1);
ds_NTFnum = ds_NTF.k*poly(cell2mat(ds_NTF.z));
ds_NTFden = poly(cell2mat(ds_NTF.p));

% All-Pole STF
ds_STFnum = [1];
ds_STFnum = ds_STFnum/sum(ds_STFnum)*sum(ds_NTFden);

% NTF Frequency Analysis
[ds_NTFjw,w] = freqz(ds_NTFnum,ds_NTFden,N);
[ds_STFjw,w] = freqz(ds_STFnum,ds_NTFden,N);
ds_NTFmag    = db( abs(ds_NTFjw) );
ds_STFmag    = db( abs(ds_STFjw) );

% In-Band Peak and RMS Gain
ds_NTFmag_stopband     = ds_NTFmag(mod(order,2)+1:fB);
ds_NTFmag_stopband_avg = mean(ds_NTFmag_stopband);
ds_NTFmag_stopband_max =  max(ds_NTFmag_stopband);

% Assign Dynamic Range to maximum DelSig In-Band Gain Value
DR = 1.1*abs(ds_NTFmag_stopband_max);

%% HTGA NTF

myNTF = HTGA(@HP_obj,@HP_init);

% Evolved NTF
theta = myNTF.F_s.best_chr;
[myNTFnum,myNTFden] = theta2poly(theta);

% [myNTFnum,myNTFden] = cheby2(order,DR,1/OSR,'high'); 

% Calculate Evolved NTF/STF Frequency Response
[myNTFjw,w] = freqz(myNTFnum,myNTFden,N);
myNTFmag    = db( abs(myNTFjw) );

%% Chebyshev II NTF (baseline)

[b,a]   = cheby2(order,DR,1/OSR,'high'); 

% All-Pole STF
STFnum = [ 1 ];
STFnum = STFnum/sum(STFnum)*sum(a);

% Calculate NTF/STF Frequency Responses
[NTFjw,w] = freqz(b,a,N);
[STFjw,w] = freqz(STFnum,a,N);
NTFmag    = db( abs(NTFjw) );
STFmag    = db( abs(STFjw) );

%% Learning Curve Plot

best_cost  = myNTF.metric_s.best_cost;
avg_cost   = myNTF.metric_s.avg_cost;
worst_cost = myNTF.metric_s.worst_cost;

figure(1)
clf
n = 1:length(best_cost);
subplot(211)
semilogy(n,worst_cost,n,avg_cost)
legend('Worst Cost','Average Cost','Location','northeast')
title('HTGA Population Fitness')
ylabel('J_{\it{x}}')
subplot(212)
plot(n,best_cost)
legend('Best Cost','Location','northeast')
xlabel('Generations (i)')
ylabel('J_{\it{x}}')

%% NTF Magnitude Response Plots

figure(2)
clf
plot(w/pi,NTFmag,'--b')
hold on
plot(w/pi,myNTFmag,'b','LineWidth',2);
plot(w/pi,ds_NTFmag,'-.r');
axis([0 0.5 -110 10]);
grid on
title('\Delta\SigmaM NTF and STF Magnitude Response (dB)')
xlabel('Normalized Frequency (x\pi rad/sample)')
ylabel('Magnitude (dB)')
legend('Traditional NTF','Evolved NTF','Delsig NTF','Location','southeast')

% Cut-Off Frequency Line
line([ 1/OSR 1/OSR ],[ -DR 0 ],'LineStyle','--','LineWidth',2,'Color','m');
text(1/OSR,0,' \leftarrow OSR^{-1}','FontSize',14); 

% Dynamic Range Line
line([ 0 1/OSR ],[ -DR -DR ],'LineStyle','--','LineWidth',2,'Color','m');
text(1/OSR,-DR,' \leftarrow DR','FontSize',14);
 
%% Evolved NTF Pole-Zero Diagram

figure(3)
clf

% HTGA PZ Plot
zplane(myNTFnum,myNTFden)
title('HTGA \Delta\SigmaM NTF Pole-Zero Diagram')

N_old = N;

% disp('NTF Done ... press any key to continue.')
% pause

%% DelSig ToolBox Analysis
%  The following analyses are done utilizing the routines available in the
% DelSig Toolbox.
% * Discrete-Time Simulation
% * Hann Windowed Output Spectrum
% * Predicted & Simulated SNR vs. Input Magnitude

myNTFzpk = zpk(roots(myNTFnum),roots(myNTFden),sum(myNTFnum),-1)
N = 8192; fB = ceil(N/(2*OSR)); f = 85;
u = 0.5*sin(2*pi*f/N*[0:N-1]);
v = simulateDSM(u,myNTFzpk);

[snr_pred,amp] = predictSNR(myNTFzpk,64);
[snr,amp] = simulateSNR(myNTFzpk,64);

figure(4)
clf
t = 0:250;
stairs(t, u(t+1),'r')
hold on
stairs(t, v(t+1))
axis([0 max(t) -1.2 1.2])
ylabel('u, v');

figure(5)
clf
spec = fft(v.*hann(N))/(N/4);
plot(linspace(0,1,N/2), dbv(spec(1:N/2)));
ylabel('dB');
snr = calculateSNR(spec(1:fB),f);
s = sprintf('SNR = %4.1fdB\n',snr);
text(0.25,-90,s);

figure(6)
clf
plot(amp,snr_pred,'b',amp,snr,'gs')
figureMagic([-100 0],10,1,[0 100], 10, 1);
xlabel('Input Level, dB')
title('SNR curve');
s = sprintf('peak SNR = %4.1fdB\n', max(snr));
text(-49,15,s)

disp('DelSig Analysis Done ... press any key to continue.')
pause

%% HTGA STF
N = N_old;

if (makeSTF == 1)
    
    % HTGA Evolved STF
    mySTF = HTGA(@LPSTF_obj,@LPSTF_init);
    gamma = mySTF.F_s.best_chr;
    mySTFnum = gamma2poly(gamma);
    mySTFnum = mySTFnum/sum(mySTFnum)*sum(myNTFden);

else
    
    % All-Pole STF
    mySTFnum = [ 1 ];
    mySTFnum = mySTFnum/sum(mySTFnum)*sum(myNTFden);

end

% Calculate Evolved STF Frequency Response
[mySTFjw,w] = freqz(mySTFnum,myNTFden,N);
mySTFmag    = db( abs(mySTFjw) );

figure(7);clf
plot(w/pi,myNTFmag,w/pi,mySTFmag)
grid on
title('HTGA Evolved \Delta\SigmaM Magnitude Response (dB)')
xlabel('Normalized Frequency (x \Pi rad/sample)')
ylabel('Magnitude (dB)')
legend('NTF','STF','Location','southeast')
axis([0 1 -2*DR 10])

figure(8);clf

subplot(121)
zplane(mySTFnum,myNTFden)
axis([-1.2 1.2 -1.2 1.2])
title('HTGA STF Pole-Zero Diagram')

subplot(122)
zplane(myNTFnum,myNTFden)
axis([-1.2 1.2 -1.2 1.2])
title('HTGA NTF Pole-Zero Diagram')