function [NTF,STF,DR] = delsigLP(order,OSR)

%% Initialization

N = 2^10; % 1024 point FFT

%% DelSig Toolbox NTF
ds_NTF    = synthesizeNTF(order,OSR,1);
ds_NTFnum = ds_NTF.k*poly(cell2mat(ds_NTF.z));
ds_NTFden = poly(cell2mat(ds_NTF.p));

%% All-Pole STF

ds_STFnum = [1 zeros(1,order)];
ds_STFnum = ds_STFnum/sum(ds_STFnum)*sum(ds_NTFden);

%% NTF Frequency Analysis
[ds_NTFjw,w] = freqz(ds_NTFnum,ds_NTFden,N);
ds_NTFmag    = db( abs(ds_NTFjw) );

%% In-Band Peak and RMS Gain

fB    = ceil(N/(2*OSR));
ds_NTFmag_stopband     = ds_NTFmag(mod(order,2)+1:fB);
ds_NTFmag_stopband_max =  max(ds_NTFmag_stopband);

% Assign Dynamic Range to maximum DelSig In-Band Gain Value
DR = abs(ds_NTFmag_stopband_max);

NTF.num = ds_NTFnum;
NTF.den = ds_NTFden;
STF.num = ds_STFnum;
STF.den = ds_NTFden;
