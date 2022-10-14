set(0,'defaulttextinterpreter','none')

close all
clear all

ORDER =    1;
OSR   =   8;
DR    =   10;
N     = 2^14;
F_C  = 1/OSR;

[NTF_1,STF_1] = chebyLP(1,OSR,DR);
[NTF_2,STF_2] = chebyLP(2,OSR,DR);

% NTF_2.num = [0 0 1];
% NTF_2.den = [1 -2 1];

[NTF_1_jw,w] = freqz(NTF_1.num,NTF_1.den,N);
% [STF_1_jw,w] = freqz(STF_1.num,STF_1.den,N);

[NTF_2_jw,w] = freqz(NTF_2.num,NTF_2.den,N);
% [STF_2_jw,w] = freqz(STF_2.num,STF_2.den,N);

NTF_1_dB = db(abs(NTF_1_jw));
% STF_1_dB = db(abs(STF_1_jw));

NTF_2_dB = db(abs(NTF_2_jw));
% STF_2_dB = db(abs(STF_2_jw));

% NTF Magnitude Repsonsess
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
% hold on
plot(w./pi,NTF_1_dB,'LineWidth',1.5)
hold on
plot(w./pi,NTF_2_dB,'LineStyle','--','Color','r')
% plot(w./pi,AA_dB,'k')
% title('Discrete-Time Low-Pass $\Delta\Sigma$ Modulator Frequency Response')
xlabel('Frequency (rad/sample)')
ylabel('Magnitude (dB)')
axis([0 0.5 -(DR+20) 10])
% set(gca,'XTick',0:1/4:1/2)
% set(gca,'XTickLabel',{'0','$\pi/4$','$\pi/2$'})

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

figure(2);
zplane(NTF_2.num,NTF_2.den)



