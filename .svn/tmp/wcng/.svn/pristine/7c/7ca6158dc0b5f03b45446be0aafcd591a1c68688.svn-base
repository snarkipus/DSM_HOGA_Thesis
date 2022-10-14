function [hogaNTF,hogaSTF] = growLP(config_s)

global order

hogaOutput = HOGA(config_s);
[hogaNTF.num,hogaNTF.den] = theta2poly(hogaOutput.F_s.best_chr);

% All-Pole STF
STFnum = [ 1 zeros(1,order)];
STFnum = STFnum/sum(STFnum)*sum(hogaNTF.den);

hogaSTF.num = STFnum;
hogaSTF.den = hogaNTF.den;

