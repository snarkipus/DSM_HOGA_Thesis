function [NTF,STF] = chebyLP(order,OSR,DR)

%% Chebyshev II NTF (baseline)

[b,a]   = cheby2(order,DR,1/OSR,'high'); 

% All-Pole STF
STFnum = [ 1 zeros(1,order)];
STFnum = STFnum/sum(STFnum)*sum(a);

NTF.num = b;
NTF.den = a;
STF.num = STFnum;
STF.den = NTF.den;

[NTF_num,den]    = cheby2(Order,R,2*pi*1.4e6,'high','s');
[STF_num,STFden] = cheby2(Order,90,2*pi*3e7,'s');
STF_num          = STF_num/STF_num(end)*den(end);