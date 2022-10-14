clear all
close all

%% Initialization

global OSR DR order

tic

% DSM Information
order   =    4;
OSR     =   32;
FFT_LEN =  2^10;

% GA Coefficients and Settings
config_s.P               =   100; % Population Size (# of chromosomes)
config_s.CR              =   0.2; % Cross-Over Rate
config_s.MR              =  0.02; % Mutation Rate
config_s.max_iter        = 10000; % Maximum Number of Iterations
config_s.min_iter        =  1000; % Minimum Number of Iterations
config_s.converge_min    =    50; % Minimum Number of Iterations for Convergence
config_s.epsilon         = 1e-20; % Misadjustment threshold

% Simulation Information
SIM_CONFIG.db_down =       3;
SIM_CONFIG.runs    =       2;
SIM_CONFIG.FFT_LEN = FFT_LEN;
SIM_CONFIG.f0      =     2/9;

% Command-Line Output
timestamp = clock;
disp(' ')
disp(['START TIME >> ' num2str(timestamp(4)) ':' num2str(timestamp(5))])
disp(' ')
 
%% DelSig DSM System Design

model_name = ['LP_' num2str(order) '_OSR_' num2str(OSR) '_ds.mat'];
file_name  = ['./RESULTS/' model_name];

try
    ds_model = open(file_name);
    disp(' ')
    disp('Loading DelSig Model ...')
    disp(' ')
    dsNTF    = ds_model.dsNTF;
    dsSTF    = ds_model.dsSTF;
    DR       = ds_model.DR + 1;
catch
    disp(' ')
    disp('Generating DelSig Model ...')
    disp(' ')

    [dsNTF,dsSTF,DR] = delsigLP(order,OSR);
    
    % Chebyshev needs to pwn Schrier
    DR = DR + 1;

    save(file_name,'dsNTF','dsSTF','DR');
end

[dsNTF.jw,w] = freqz(dsNTF.num,dsNTF.den,FFT_LEN);
[dsSTF.jw,w] = freqz(dsSTF.num,dsSTF.den,FFT_LEN);
dsNTF.mag    = db( abs(dsNTF.jw) );
dsSTF.mag    = db( abs(dsSTF.jw) );

%% Traditional DSM System Design

model_name = ['LP_' num2str(order) '_OSR_' num2str(OSR) '_ch.mat'];
file_name  = ['./RESULTS/' model_name];

try
    ch_model = open(file_name);
    disp(' ')
    disp('Loading ChebyII Model ...')
    disp(' ')
    chNTF    = ch_model.chNTF;
    chSTF    = ch_model.chSTF;
catch
    disp(' ')
    disp('Generating ChebyII Model ...')
    disp(' ')

    [chNTF,chSTF] = chebyLP(order,OSR,DR);

    save(file_name,'chNTF','chSTF');
end

% We have to a little bit better
DR = DR + 3;

[chNTF.jw,w] = freqz(chNTF.num,chNTF.den,FFT_LEN);
[chSTF.jw,w] = freqz(chSTF.num,chSTF.den,FFT_LEN);
chNTF.mag    = db( abs(chNTF.jw) );
chSTF.mag    = db( abs(chSTF.jw) );

%% DSM Simulation

sim_name  = ['LP_' num2str(order) '_OSR_' num2str(OSR) '_sim.mat'];
file_name  = ['./RESULTS/' sim_name];

disp(' ')
clobber = input('XXX Destroy previous model simulations (if exists) XXX [n]: ','s');
if isempty(clobber)
    clobber = 'n';
end

if clobber == 'n'
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
else
    disp(' ')
    disp('Simulating Models ...')
    disp(' ')

    dsOutput = DTsimulation(dsNTF,dsSTF,OSR,SIM_CONFIG);
    chOutput = DTsimulation(chNTF,chSTF,OSR,SIM_CONFIG);
    save(file_name,'dsOutput','chOutput');
end

%% HTGA DSM System Design

HOGA_done = 'FALSE';

model_name = ['LP_' num2str(order) '_OSR_' num2str(OSR) '_hoga.mat'];
file_name  = ['./RESULTS/' model_name];

% //////////////////////
% /// DSM EVOLUTION LOOP

while ( strcmp(HOGA_done,'FALSE') )

    % //////////////////////
    % /// NTF EVOLUTION LOOP

    disp(' ')
    clobber = input('XXX Destroy previous NTF model (if exists) XXX [n]: ','s');
    if isempty(clobber)
        clobber = 'n';
    end

    NTF_done  = 'FALSE';

    while ( strcmp(NTF_done,'FALSE') )

        % PRESERVE PREVIOUS MODEL (IF AVAILABLE)
        if clobber == 'n'
            try
                hoga_model = open(file_name);
                hogaNTF    = hoga_model.hogaNTF;
                hogaSTF    = hoga_model.hogaSTF;
            catch
                disp(' ')
                disp('Evolving NTF ...')
                disp(' ')

                config_s.init_fcn = @HP_init;
                [hogaNTF,hogaSTF] = growLP(config_s);

                save(file_name,'hogaNTF','hogaSTF');
            end

        else
            % CREATE NEW MODEL (DESTRUCTIVE)
            disp(' ')
            disp('Evolving NTF ...')
            disp(' ')

            config_s.init_fcn = @HP_init;
            [hogaNTF,hogaSTF] = growLP(config_s);

            save(file_name,'hogaNTF','hogaSTF');
        end


        [hogaNTF.jw,w] = freqz(hogaNTF.num,hogaNTF.den,FFT_LEN);
        [hogaSTF.jw,w] = freqz(hogaSTF.num,hogaSTF.den,FFT_LEN);
        hogaNTF.mag    = db( abs(hogaNTF.jw) );
        hogaSTF.mag    = db( abs(hogaSTF.jw) );

        % //////////////////////
        % /// VisualizeNTF

        % Magnitude Response Comparison
        figure(1);clf
        plot(w/pi,chNTF.mag,'g')
        hold on
        plot(w/pi,hogaNTF.mag,'b','LineWidth',2);
        plot(w/pi,dsNTF.mag,'r');
        axis([0 0.5 -110 10]);
        grid on
        title('\Delta\SigmaM NTF Magnitude Response Comparison')
        xlabel('Normalized Frequency (\times\pi rad/sample)')
        ylabel('Magnitude (dB)')
        legend('Chebyshev','HOGA','Delsig','Location','southeast')
%         line([ 1/OSR 1/OSR ],[ -DR+3     0 ],'LineStyle','--','LineWidth',2,'Color','m');
%         line([     0 1/OSR ],[ -DR+3 -DR+3 ],'LineStyle','--','LineWidth',2,'Color','m');

        % Pole-Zero Diagram
        figure(2);clf
        zplane(hogaNTF.num,hogaNTF.den)
        title('hoga \Delta\SigmaM NTF Pole-Zero Diagram')

        disp(' ')
        loop = input('Repeat NTF Evolution[n]: ','s');
        if isempty(loop)
            loop = 'n';
        end

        if loop == 'n'
            NTF_done = 'TRUE';
        else
            clobber = 'y';
        end

    end

    % /// END EVOLUTION LOOP
    % //////////////////////

    disp(' ')
    clobber = input('XXX Destroy previous STF model (if exists) XXX [n]: ','s');
    if isempty(clobber)
        clobber = 'n';
    end

    % //////////////////////
    % /// STF EVOLUTION LOOP

    STF_done  = 'FALSE';

    while ( strcmp(STF_done,'FALSE') )

        % CREATE NEW MODEL (DESTRUCTIVE)
        if clobber ~= 'n'
            disp(' ')
            disp('Evolving STF ...')
            disp(' ')

            config_s.init_fcn = @LPSTF_init;
            hogaSTF = growLP_STF(config_s,hogaNTF.den);

            save(file_name,'hogaNTF','hogaSTF');
        end


        [hogaNTF.jw,w] = freqz(hogaNTF.num,hogaNTF.den,FFT_LEN);
        [hogaSTF.jw,w] = freqz(hogaSTF.num,hogaSTF.den,FFT_LEN);
        hogaNTF.mag    = db( abs(hogaNTF.jw) );
        hogaSTF.mag    = db( abs(hogaSTF.jw) );

        % //////////////////////
        % /// VisualizeSTF

        % Magnitude Response Comparison
        figure(3);clf
        plot(w/pi,chSTF.mag,'g')
        hold on
        plot(w/pi,hogaSTF.mag,'b','LineWidth',2);
        plot(w/pi,dsSTF.mag,'r');
        axis([0 0.5 -110 10]);
        grid on
        title('\Delta\SigmaM STF Magnitude Response Comparison')
        xlabel('Normalized Frequency (x\pi rad/sample)')
        ylabel('Magnitude (dB)')
        legend('Classical','HOGA','Delsig','Location','southeast')

        % Pole-Zero Diagram
        figure(4);clf
        zplane(hogaSTF.num,hogaSTF.den)
        title('hoga \Delta\SigmaM STF Pole-Zero Diagram')

        disp(' ')
        loop = input('Repeat STF Evolution[n]: ','s');
        if isempty(loop)
            loop = 'n';
        end

        if loop == 'n'
            STF_done = 'TRUE';
        else
            clobber = 'y';
        end

    end

    % /// END EVOLUTION LOOP
    % //////////////////////

    % hoga Simulation
    hogaOutput = DTsimulation(hogaNTF,hogaSTF,OSR,SIM_CONFIG);

    table1 = ...
        {'SIMULATION', 'RESULTS' , ['ORDER: ' num2str(order)], ['OSR: ' num2str(OSR)];...
        '**********'  ,'*******' ,                 '********',             '********';...
        'TYPE'   ,   'SNR[dB]'   ,               'DR[dB]'    ,            'FLOOR[dB]';...
        'DelSig' ,   dsOutput.SNR,                dsOutput.DR,          dsOutput.SFDR;...
        'ChebyII',   chOutput.SNR,                chOutput.DR,          chOutput.SFDR;...
        'hoga'   , hogaOutput.SNR,              hogaOutput.DR,        hogaOutput.SFDR};

    disp(' ')
    disp(table1)

    % Output Spectrum Comparison
    figure(5); clf
    n = (1:length(dsOutput.mag)).*(4/(FFT_LEN*OSR));
    plot(n,chOutput.mag,'g')
    hold on
    plot(n,dsOutput.mag,'r')
    plot(n,hogaOutput.mag,'b','LineWidth',2)
    title({'\Delta\SigmaM Output Spectrum Comparison';...
        ['Order:' num2str(order)...
        ' OSR:' num2str(OSR)...
        ' Runs:' num2str(SIM_CONFIG.runs)...
        ' FFT:' num2str(SIM_CONFIG.FFT_LEN)]})
    axis([0 2/OSR min(dsOutput.mag) 10])
    legend('Chebyshev','DelSig','HOGA')
    ylabel('Magnitude (dB)')
    xlabel('Normalized Frequency (x\pi rad/sample)')

    figure(6); clf
    plot(w/pi,chNTF.mag,'g');
    hold on
    plot(w/pi,dsNTF.mag,'r');
    plot(w/pi,hogaNTF.mag,'b','LineWidth',2);
    plot(w/pi,chSTF.mag,'g');
    plot(w/pi,dsSTF.mag,'r');
    plot(w/pi,hogaSTF.mag,'b','LineWidth',2);
    axis([0 0.5 -110 10]);
    axis square
    grid on
    title({'\Delta\SigmaM Magnitude Response Comparison';
            ['Order:' num2str(order)...
            ' OSR:' num2str(OSR)]})
    xlabel('Normalized Frequency (\times\pi rad/sample)')
    ylabel('Magnitude (dB)')
    legend('Chebyshev','Delsig','HOGA','Location','southeast')   
    
    disp(' ')
    loop = input('Repeat hoga Evolution[n]: ','s');
    if isempty(loop)
        loop = 'n';
    end

    if loop == 'n'
        HOGA_done = 'TRUE';
    else
        clobber = 'y';
    end
end

figure(7);clf

subplot(121)
    plot(w/pi,chNTF.mag,'g');
    hold on
    plot(w/pi,dsNTF.mag,'r');
    plot(w/pi,hogaNTF.mag,'b','LineWidth',2);
    plot(w/pi,chSTF.mag,'g');
    plot(w/pi,dsSTF.mag,'r');
    plot(w/pi,hogaSTF.mag,'b','LineWidth',2);
    axis([0 0.5 -110 10]);
    axis square
    grid off
    title({'\Delta\SigmaM Magnitude Response Comparison';
            ['Order:' num2str(order)...
            ' OSR:' num2str(OSR)]})
    xlabel('Normalized Frequency (\times\pi rad/sample)')
    ylabel('Magnitude (dB)')
    legend('Chebyshev','Delsig','HOGA','Location','southeast')   

subplot(122)
    n = (1:length(dsOutput.mag)).*(4/(FFT_LEN*OSR));
    plot(n,chOutput.mag,'g')
    hold on
    plot(n,dsOutput.mag,'r')
    plot(n,hogaOutput.mag,'b','LineWidth',2)
    title({'\Delta\SigmaM Output Spectrum Comparison';...
        ['Order:' num2str(order)...
        ' OSR:' num2str(OSR)...
        ' Runs:' num2str(SIM_CONFIG.runs)...
        ' FFT:' num2str(SIM_CONFIG.FFT_LEN)]})
    axis([0 2/OSR min(dsOutput.mag) 10])
    axis square
    grid off
    legend('Chebyshev','DelSig','HOGA')
    ylabel('Magnitude (dB)')
    xlabel('Normalized Frequency (x\pi rad/sample)')
    
    toc
    