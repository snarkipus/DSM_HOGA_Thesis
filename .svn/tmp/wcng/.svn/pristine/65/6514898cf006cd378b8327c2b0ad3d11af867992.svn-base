clc
clear all
close all

%% Initialization

global OSR DR order

% DSM Information
order   =     3;
OSR     =    32;
FFT_LEN =  2^10;

% Simulation Information
SIM_CONFIG.db_down =       3;
SIM_CONFIG.runs    =       3;
SIM_CONFIG.FFT_LEN = FFT_LEN;
SIM_CONFIG.f0      =     1/5;

dir_name = ['../RESULTS/ORDER_' num2str(order) '/OSR_' num2str(OSR) '/'];
% dir_name = ['../temp/RESULTS/ORDER_' num2str(order) '/OSR_' num2str(OSR) '/'];


% Command-Line Output
timestamp = clock;
disp(' ')
disp(['START TIME >> ' num2str(timestamp(4)) ':' num2str(timestamp(5))])
disp(' ')

%% DelSig DSM System Design

model_name = ['LP_' num2str(order) '_OSR_' num2str(OSR) '_ds.mat'];
file_name  = [dir_name model_name];

try
    ds_model = open(file_name);
    dsNTF    = ds_model.dsNTF;
    dsSTF    = ds_model.dsSTF;
    DR       = ds_model.DR;
catch
    disp(' ')
    disp('Generating DelSig Model ...')
    disp(' ')

    [dsNTF,dsSTF,DR] = delsigLP(order,OSR);

    [dsNTF.jw,w] = freqz(dsNTF.num,dsNTF.den,FFT_LEN);
    [dsSTF.jw,w] = freqz(dsSTF.num,dsSTF.den,FFT_LEN);
    dsNTF.mag    = db( abs(dsNTF.jw) );
    dsSTF.mag    = db( abs(dsSTF.jw) );

    save(file_name,'dsNTF','dsSTF','DR');
end

[dsNTF.jw,w] = freqz(dsNTF.num,dsNTF.den,FFT_LEN);
dsNTF.mag    = db( abs(dsNTF.jw) );

%% Traditional DSM System Design

model_name = ['LP_' num2str(order) '_OSR_' num2str(OSR) '_ch.mat'];
file_name  = [dir_name model_name];

try
    ch_model = open(file_name);
    chNTF    = ch_model.chNTF;
    chSTF    = ch_model.chSTF;
catch
    disp(' ')
    disp('Generating ChebyII Model ...')
    disp(' ')
    
    DRch = DR + 2;    
    
    [chNTF,chSTF] = chebyLP(order,OSR,DRch);

    [chNTF.jw,w] = freqz(chNTF.num,chNTF.den,FFT_LEN);
    [chSTF.jw,w] = freqz(chSTF.num,chSTF.den,FFT_LEN);
    chNTF.mag    = db( abs(chNTF.jw) );
    chSTF.mag    = db( abs(chSTF.jw) );

    save(file_name,'chNTF','chSTF');
end

%% DSM Simulation

sim_name  = ['LP_' num2str(order) '_OSR_' num2str(OSR) '_sim.mat'];
file_name = [dir_name sim_name];

try
    Output   = open(file_name);
    dsOutput = Output.dsOutput;
    chOutput = Output.chOutput;
catch
    disp(' ')
    disp('Simulating Models ...')
    disp(' ')

    dsOutput = DTsimulation(dsNTF,dsSTF,OSR,SIM_CONFIG);
    chOutput = DTsimulation(chNTF,chSTF,OSR,SIM_CONFIG);
    save(file_name,'dsOutput','chOutput');
end

%% HTGA DSM System Design

DR = DR + 3;

HTGA_done = 'FALSE';

model_name = ['LP_' num2str(order) '_OSR_' num2str(OSR) '_htga.mat'];
file_name  = [dir_name model_name];

% //////////////////////
% /// DSM EVOLUTION LOOP

while ( strcmp(HTGA_done,'FALSE') )

    disp(' ')
    disp('Evolving NTF ...')
    disp(' ')

    [htgaNTF,htgaSTF] = growLP();
    [htgaNTF.jw,w] = freqz(htgaNTF.num,htgaNTF.den,FFT_LEN);
    htgaNTF.mag    = db( abs(htgaNTF.jw) );

%     disp(' ')
%     disp('Evolving STF ...')
%     disp(' ')
% 
%     htgaSTF = growLP_STF(htgaNTF.den);
    [htgaSTF.jw,w] = freqz(htgaSTF.num,htgaSTF.den,FFT_LEN);
    htgaSTF.mag    = db( abs(htgaSTF.jw) );

    save(file_name,'htgaNTF','htgaSTF');

    % HTGA Simulation
    htgaOutput = DTsimulation(htgaNTF,htgaSTF,OSR,SIM_CONFIG);

    table1 = ...
        {'SIMULATION', 'RESULTS' , ['ORDER: ' num2str(order)], ['OSR: ' num2str(OSR)];...
        '**********'  ,'*******' ,                 '********',             '********';...
        'TYPE'   ,   'SNR[dB]'   ,               'DR[dB]'    ,            'FLOOR[dB]';...
        'DelSig' ,   dsOutput.SNR,                dsOutput.DR,          dsOutput.SFDR;...
        'ChebyII',   chOutput.SNR,                chOutput.DR,          chOutput.SFDR;...
        'HTGA'   , htgaOutput.SNR,              htgaOutput.DR,        htgaOutput.SFDR};

    disp(' ')
    disp(table1)

    % Output Spectrum Comparison
    figure(2); clf
    n = (1:length(dsOutput.mag)).*(4/(FFT_LEN*OSR));
    plot(n,dsOutput.mag,n,chOutput.mag,n,htgaOutput.mag)
    title({'\Delta\SigmaM Output Spectrum Comparison';...
        ['Order:' num2str(order)...
        ' OSR:' num2str(OSR)...
        ' Runs:' num2str(SIM_CONFIG.runs)...
        ' FFT:' num2str(SIM_CONFIG.FFT_LEN)]})
    axis([0 2/OSR min(dsOutput.mag) 10])
    legend('DelSig','ChebyII','HTGA')
    ylabel('Magnitude (dB)')
    xlabel('Normalized Frequency (x\pi rad/sample)')

    if (OSR == 128 && isreal(htgaOutput.SNR) && htgaOutput.SNR > 0)
        HTGA_done = 'TRUE';
    elseif ( htgaOutput.SNR > dsOutput.SNR && htgaOutput.DR > dsOutput.DR )
        HTGA_done = 'TRUE';
    end
    
end

figure(1);clf

% //////////////////////
% /// VisualizeNTF

% Magnitude Response Comparison
subplot(221)
plot(w/pi,chNTF.mag,'--b')
hold on
plot(w/pi,htgaNTF.mag,'b','LineWidth',2);
plot(w/pi,dsNTF.mag,'-.r');
axis([0 0.5 -110 10]);
grid on
title('\Delta\SigmaM NTF Magnitude Response Comparison')
xlabel('Normalized Frequency (x\pi rad/sample)')
ylabel('Magnitude (dB)')
legend('Traditional NTF','Evolved NTF','Delsig NTF','Location','southeast')
line([ 1/OSR 1/OSR ],[ -DR+3     0 ],'LineStyle','--','LineWidth',2,'Color','m');
line([     0 1/OSR ],[ -DR+3 -DR+3 ],'LineStyle','--','LineWidth',2,'Color','m');

% Pole-Zero Diagram
subplot(222)
zplane(htgaNTF.num,htgaNTF.den)
title('HTGA \Delta\SigmaM NTF Pole-Zero Diagram')

% //////////////////////
% /// VisualizeSTF

% Magnitude Response Comparison
subplot(223)
plot(w/pi,chSTF.mag,'--b')
hold on
plot(w/pi,htgaSTF.mag,'b','LineWidth',2);
plot(w/pi,dsSTF.mag,'-.r');
axis([0 0.5 -110 10]);
grid on
title('\Delta\SigmaM STF Magnitude Response Comparison')
xlabel('Normalized Frequency (x\pi rad/sample)')
ylabel('Magnitude (dB)')
legend('Traditional STF','Evolved STF','Delsig STF','Location','southeast')

% Pole-Zero Diagram
subplot(224)
zplane(htgaSTF.num,htgaSTF.den)
title('HTGA \Delta\SigmaM STF Pole-Zero Diagram')
