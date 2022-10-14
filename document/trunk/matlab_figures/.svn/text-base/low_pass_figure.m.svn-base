set(0,'defaulttextinterpreter','tex')

close all
clear all

ORDER =    5;
OSR   =   32;
DR    =   70;
N     = 2^14;
F_C  = 1/OSR;

% Optimal Equiripple FIR Design
pb_ripple = 0.1;
sb_ripple = 80;
f = [0.8*F_C 1.25*F_C];
a = [1 0];
dev = [(10^(pb_ripple/20)-1)/(10^(pb_ripple/20)+1)  10^(-sb_ripple/20)];
[FIR_ORDER,f0,a0,w]=firpmord(f,a,dev);
b = firpm(FIR_ORDER,f0,a0,w);

[AA_jw,w] = freqz(b,[1],N);
AA_dB = db(abs(AA_jw));

[NTF,STF] = chebyLP(ORDER,OSR,DR);

[NTF_jw,w] = freqz(NTF.num,NTF.den,N);
[STF_jw,w] = freqz(STF.num,STF.den,N);

NTF_dB = db(abs(NTF_jw));
STF_dB = db(abs(STF_jw));

% NTF Magnitude Repsonses
figure(1);clf

FUDGE_X = .05*(1/OSR);

% % Operational Region
% rectangle('Position',[0+FUDGE_X,-DR+0.5,1/OSR-2*FUDGE_X,DR-1],...
%           'Curvature',[0.05],...
%           'LineWidth',1.5)
% text(1/(2*OSR),-DR/2,'$\Delta\Sigma$M OPERATIONAL REGION',...
%     'rotation',90,...
%     'HorizontalAlignment','center')
% 
% hold on

% % NTF Annotation
% text(4/10,3,'Noise Shaping Reponse: $\mathcal{T}_n[\cdot]$',...
%     'rotation',0,...
%     'HorizontalAlignment','center')
% 
% % % STF Annotation
% text(4/10,-60,'Signal Shaping Reponse: $\mathcal{T}_s[\cdot]$',...
%     'rotation',-20,...
%     'HorizontalAlignment','center')
hold on
plot(w./pi,NTF_dB,'LineWidth',1)
plot(w./pi,STF_dB,'LineWidth',1,'Color','r','LineStyle','--')
plot(w./pi,AA_dB,'LineWidth',1,'Color','k')
% title('Discrete-Time Low-Pass $\Delta\Sigma$ Modulator Frequency Response')
xlabel('Normalized Frequency (\times \pi rad/sample)')
ylabel('Magnitude (dB)')
axis([0 0.5 -(DR+20) 5])
% set(gca,'XTick',0:1/4:1/2)
% set(gca,'XTickLabel',{'0','$\pi/4$','$\pi/2$'})
legend('NTF(z)','STF(z)','H_d(z)')
% legend('$\textrm{NTF}(\omega)$','$\textrm{STF}(\omega)$','$\textrm{LP}(\o
% mega)$')
% legend('Position',[0.5701 0.5036 0.3312 0.1577])
% legend('boxoff')
% legend('Position',[ ])


% % Cut-Off Frequency Line
% line([ 1/OSR 1/OSR ],[ -DR 0 ],'LineStyle','--','LineWidth',2,'Color','m');
% text(1/OSR,0,' $\leftarrow \textrm{OSR}^{-1}$','FontSize',14); 
% 
% % Dynamic Range Line
% line([ 0 1/OSR ],[ -DR -DR ],'LineStyle','--','LineWidth',2,'Color','m');
% text(1/OSR,-DR,' $\leftarrow \textrm{DR}$','FontSize',14);
