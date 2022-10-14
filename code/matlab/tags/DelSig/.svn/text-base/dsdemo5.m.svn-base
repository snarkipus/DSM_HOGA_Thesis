% Demonstrate the simulateESL function
echo off;
if exist('LiveDemo','var') == 0
    LiveDemo=0;
end
clc;
fprintf(1,'\t\tMismatch-Shaping Unit-Element DAC\n\n');

% Specify the modulator NTF, the mismatch-shaping TF, and the number of elements
echo on;
ntf = synthesizeNTF(3,[],[],4);
M = 16;
sigma_d = 0.01;		% 1% mismatch
mtf1 = zpk(1,0,1,1);	%First-order shaping
echo off;	
A = 1/sqrt(2);		% Test tone amplitude, relative to full-scale.
f = 0.33;			% Test tone frequency, relative to fB. 
					% (Will be adjusted to be an fft bin)

if LiveDemo
    N = 2^12;
else
    N = 2^14;
end
fin = round(f*N/(2*12));
w = (2*pi/N)*fin;
echo on;
u = M*A*sin(w*[0:N-1]);
v = simulateDSM(u,ntf,M+1);	% M unit elements requires an M+1-level quant.
v = (v+M)/2;				% scale v to [0,M]
sv1 = simulateESL(v,mtf1,M);
echo off

figure(1); clf
T = 20;
subplot(211);
plotUsage(thermometer(v(1:T),M));
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','Element Usage');
title('Thermometer-Coding')
subplot(212);
plotUsage(sv1(:,1:T));
title('First-Order Shaping');
if LiveDemo
    set(1,'position',[9 204 330 525]);
    changeFig(18,.5,1);
    pause
end

ideal = v;

% DAC element values
e_d = randn(M,1); 
e_d = e_d - mean(e_d);
e_d = sigma_d * e_d/std(e_d);
ue = 1 + e_d;

% Convert v to analog form, assuming no shaping
thermom = zeros(M+1,1);
for i=1:M
    thermom(i+1) = thermom(i) + ue(i);
end
conventional = thermom(v+1)';

% Convert sv to analog form 
dv1 = ue' * sv1;

window = hann(N);
spec = fft(ideal.*window)/(M*N/8);
spec0 = fft(conventional.*window)/(M*N/8);
spec1 = fft(dv1.*window)/(M*N/8);

figure(2); clf
plotSpectrum(spec0,fin,'r');
hold on;
plotSpectrum(spec1,fin,'b');
plotSpectrum(spec,fin,'g');
axis([1e-3 0.5 -200 0]);
x1 = 2e-3; x2=1e-2; y0=-180; dy=dbv(x2/x1); y3=y0+3*dy;
plot([x1 x2 x2 x1],[y0 y0 y3 y0],'k')
text(x2, (y0+y3)/2,' 60 dB/decade')
hold off;
grid;
ylabel('PSD');
xlabel('Normalized Frequency');
if LiveDemo
    figure(1);
    set(1,'position',[9 427 200 300]);
    changeFig;
    subplot(211)
    xlabel('')
    figure(2);
    set(2,'position',[237 288 553 439]);
    changeFig(18,2,12);
    legend('thermom.','1^{st}-Order','Ideal');
    set(gcf,'NumberTitle','off'); 
    set(gcf,'Name','Output Spectra');
	pause
    changeFig;
    legend;
    set(2,'position',[238 427 325 300]);
end
legend('thermometer','rotation','ideal DAC');
fprintf(1,'Paused.\n');
pause

%Now repeat the above for second-order shaping
echo on;
mtf2 = zpk([ 1 1 ], [ 0 0 ], 1, 1);	%Second-order shaping
sv2 = simulateESL(v,mtf2,M);
echo off;

figure(1); clf
T = 20;
subplot(211);
plotUsage(sv1(:,1:T));
title('First-Order Shaping');
subplot(212);
plotUsage(sv2(:,1:T));
title('Second-Order Shaping');
if LiveDemo
    set(1,'position',[9 204 330 525]);
    changeFig(18,.5,1);
    pause
end

dv2 = ue' * sv2;
spec2 = fft(dv2.*window)/(M*N/8);

figure(2); clf
plotSpectrum(spec1,fin,'r');
hold on;
plotSpectrum(spec2,fin,'b');
plotSpectrum(spec,fin,'g');
axis([1e-3 0.5 -200 0]);
x1 = 2e-3; x2=1e-2; y0=-180; dy=dbv(x2/x1); y3=y0+3*dy;
plot([x1 x2 x2 x1],[y0 y0 y3 y0],'k')
text(x2, (y0+y3)/2,' 60 dB/decade')
legend('1^{st}-Order','2^{nd}-Order','Ideal');
hold off;
grid;
xlabel('Normalized Frequency');
ylabel('PSD');
if LiveDemo
    figure(1);
    set(1,'position',[9 427 200 300]);
    changeFig;
    subplot(211)
    xlabel('')
    figure(2);
    set(2,'position',[237 288 553 439]);
    changeFig(18,2,12);
    legend('1^{st}-Order','2^{nd}-Order','Ideal');
	pause
    changeFig;
    legend('1^{st}-Order','2^{nd}-Order','Ideal');
    set(2,'position',[238 427 325 300]);
end

