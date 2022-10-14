function mod = dsexample2
% Design example for an 8th-order binary bandpass modulator.

% Altogether too much of the code is greared toward making the graphs
% look 'pretty'.

format compact;
J = 1i;

clc
fprintf(1,'\t\t\t8th-Order Bandpass Example\n\n');

% Design parameters
order = 8;
R = 32;
nlev = 2;
f0 = 0.125;
Hinf = 1.5;

% NTF Synthesis
fprintf(1,'Doing NTF synthesis... ');
H = synthesizeNTF(order,R,2,Hinf,f0);		%Optimized zero placement
fprintf(1,'Done.\n');

figure(1); clf
set(1,'name','Poles and Zeros');
set(1,'numbertitle','off');
set(1,'position',[0 359 200 200]);
plotPZ(H,{'r','g'});
axis([-1 1 -1 1]); axis('square');
set(gca,'position', [0.1 0.07, 0.9, 0.85]);
drawnow;

figure(2); clf
set(2,'name','NTF Magnitude Response')
set(2,'numbertitle','off');
set(2,'position',[224 359 300 200]);
f = [linspace(0,f0-0.375/R,50) linspace(f0-0.375/R,f0+0.375/R,100) ...
     linspace(f0+.375/R,0.5,50)];
z = exp(J*2*pi*f);
magH = dbv(evalTF(H,z));
NTFMagHandle = subplot('position', [0.1 0.58, 0.87, 0.4]);
plot(f,magH,'r');
hold on;
plot([f0-1/(4*R) f0+1/(4*R)], -80*[1 1],'-g','linewidth',3);
axis([0 0.5 -80 10]);
grid on;
text(0.25,-70,'Normalized frequency (1\rightarrow f_s)','hor','center');

f = linspace(-0.75,0.75,300); z = exp(J*2*pi*(f/(2*R)+f0));
magH = dbv(evalTF(H,z));
subplot('position', [0.1 0.07, 0.87, 0.4]);
plot(f,magH,'-r');
axis([-.75 .75 -80 -30]);
grid on
sigma_H = dbv(rmsGain(H,f0-1/(4*R),f0+1/(4*R)));
hold on;
plot([-0.5 0.5], sigma_H*[1 1]);
plot([-0.5 0.5], sigma_H*[1 1],'o');
plot([-0.5 0.5], -80*[1 1],'-g','linewidth',3);
hold off;
text( 0, sigma_H+5, sprintf('in-band RMS gain = %5.0fdB',sigma_H),'hor','center');
text(0,-74,'Normalized frequency (1\rightarrow f_B)','hor','center');
drawnow;

% SNR Simulation (stability verification)
fprintf(1,'Doing SNR simulations... ');
figure(3); clf;
set(3,'name','SNR curve');
set(3,'numbertitle','off');
set(3,'position',[544 359 235 200]);
if nlev==2
    [snr_pred,amp_pred] = predictSNR(H,R,[],f0);
    plot(amp_pred,snr_pred,'-');
    hold on;
end
[snr,amp] = simulateSNR(H,R,[],f0,nlev);
fprintf(1,'Done.\n');
plot(amp,snr,'og');
grid on;
figureMagic([-80 0], 10, 2, [0 80], 10, 2);
xlabel('Input Level (dBFS)');
ylabel('SNR (dB)');
[peak_snr,peak_amp] = peakSNR(snr,amp);
msg = sprintf('OSR=%d',R);
text(-50,75,msg,'hor','center', 'vertical','middle');
msg = sprintf('peak SNR = %4.1fdB  \n@ amp=%4.1fdB  ',peak_snr,peak_amp);
text(0,15,msg,'hor','right');
set(gca,'position', [0.15 0.1, 0.82, 0.87]);
drawnow;

% Realization and dynamic range scaling
fprintf(1,'Doing dynamic range scaling... ');
form = 'CRFB';
[a,g,b,c] = realizeNTF(H,form);
b = [b(1) zeros(1,length(b)-1)];	% Use a single feed-in for the input
ABCD0 = stuffABCD(a,g,b,c,form);
[junk G] = calculateTF(ABCD0);
figure(2);
subplot(NTFMagHandle);	hold on;
f = linspace(0,0.5);	z = exp(2i*pi*f);
plot(f,dbv(evalTF(G,z)),'g')
hold off; drawnow;
% [ABCD umax] = scaleABCD(ABCD0,nlev,f0,[],[],[],1e4);
[ABCD umax] = scaleABCD(ABCD0,nlev,f0);
[a,g,b,c] = mapABCD(ABCD,form);
fprintf(1,'Done.\n');

fprintf(1,'Verifying dynamic range scaling... ');
u = linspace(0,0.95*umax,30);
N = 1e4; 
tone = (nlev-1)*cos(2*pi*(f0+1/(7*R))*[1:N]);
maxima = zeros(order,length(u));
for i = 1:length(u)
    ui = u(i);
    [v,xn,xmax] = simulateDSM( ui*tone, ABCD, nlev );
    maxima(:,i) = xmax(:);
    if any(xmax>1e2) 
	fprintf(1,'Warning, umax from scaleABCD was too high.\n');
	umax = ui;
	u = u(1:i);
	maxima = maxima(:,1:i);
    	break;
    end
end
figure(4); clf;
set(4,'name','State Maxima');
set(4,'numbertitle','off');
set(4,'position',[20 339 300 200]);
colors = hsv(order);
cmd = 'legend( handles';
handles = [];
for i = 1:order
    handles(i) = plot(u,maxima(i,:),'--','color',colors(i,:));
    if i==1
	hold on;
    end
    cmd = [ cmd ',''State ' num2str(i) '''' ];
    plot(u,maxima(i,:),'o','color',colors(i,:));
end
hold off; grid on;
text(umax/2,0.05,'input amplitude','hor','cen')
axis([ 0 umax 0 1]);
set(gca,'position', [0.1 0.07, 0.85, 0.85]);
cmd = [cmd ',4);'];
eval(cmd);
drawnow;
fprintf(1,'Done.\n');

% The next step would be to size capacitors for each integrator state based 
% on noise (kT/C) limits and area.

% Simulations modelling the effects of finite op-amp gain and capacitor 
% errors should also be performed.

mod.NTF = H;
mod.STF = G;
mod.nlev = nlev;
mod.f0 = f0;
mod.ABCD = ABCD;
mod.umax = umax;
mod.peak_snr = peak_snr;
mod.form = form;
mod.coefficients.a = a;
mod.coefficients.g = g;
mod.coefficients.b = b;
mod.coefficients.c = c;

