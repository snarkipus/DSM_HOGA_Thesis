function hogaSTF = growLP_STF(config_s,den)

global myNTFden

myNTFden = den;

hogaOutput = HOGA(config_s);

gamma = hogaOutput.F_s.best_chr;
hogaSTF.num = gamma2poly(gamma);
hogaSTF.num = hogaSTF.num/sum(hogaSTF.num)*sum(den);
    
hogaSTF.den = den;

