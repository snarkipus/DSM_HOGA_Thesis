function [htgaNTF,htgaSTF] = growLP()

global order

htgaOutput = HTGA(@HP_obj,@HP_init);
[htgaNTF.num,htgaNTF.den] = theta2poly(htgaOutput.F_s.best_chr);

% All-Pole STF
STFnum = [ 1 zeros(1,order)];
STFnum = STFnum/sum(STFnum)*sum(htgaNTF.den);

htgaSTF.num = STFnum;
htgaSTF.den = htgaNTF.den;

