%% Continous Time Signal
% % $$x(t)=A \cos(2\pi f t)$$ %
CTsig.freq   = 1000;
CTsig.per    = CTsig.freq^-1;
CTsig.amp    = 1;
CTsig.cycles = 4;
CTsig.time   = 0:CTsig.per/10000:CTsig.cycles*CTsig.per;
CTsig.signal = CTsig.amp*cos(2*pi*CTsig.freq*CTsig.time);

%% Discretized (Sampled) Signal 
% % $$x(n)=A \cos\left(2\pi f (n T_s) \right)$$ % 
DTsig.osr    = 8; 
DTsig.f_s    = DTsig.osr*2*CTsig.freq; 
DTsig.T_s    = DTsig.f_s^-1; 
DTsig.n      = 0:2*CTsig.cycles*DTsig.osr; 
DTsig.signal = CTsig.amp*cos(2*pi*CTsig.freq*(DTsig.n*DTsig.T_s)); 

%% Output 
figure(1), clf 
plot(1000*CTsig.time,CTsig.signal,'--r')
xlabel('Time [ms]')
ylabel('Amplitude')
title('Sampling')
hold on
stem(DTsig.n/(DTsig.osr*2),DTsig.signal,'fill')

%% Quantized Signal
% % $$\hat{x}(n) = \mathcal{Q}[x(n)]=x(n)+e(n)$$ %
Qsig.n      = DTsig.n;
Qsig.bits   = 2;
Qsig.L      = 2^Qsig.bits;
Qsig.scale  = DTsig.signal*(Qsig.L - 0.5);
Qsig.shift  = Qsig.scale-0.5;
Qsig.round  = round(Qsig.shift);
Qsig.signal = (Qsig.round+0.5)/Qsig.L;
Qsig.error  = DTsig.signal - Qsig.signal;

%% Zero Order Hold Output
% % $$x_{\mathrm{ZOH}}(t)=\sum_{n=-\infty}^{\infty}x(n)\frac{\sin \pi (t - n T_s)}{\pi(t - n T_s)}$$ % 
figure(2);clf
stem(DTsig.n/(DTsig.osr*2),Qsig.signal,'fill')
hold on
stairs(DTsig.n/(DTsig.osr*2),Qsig.signal,'r','LineWidth',2)
hold off
xlabel('Time [ms]')
ylabel('Amplitude')
title('Zero-Order Hold Output')
set(gca,'YGrid','on')

   