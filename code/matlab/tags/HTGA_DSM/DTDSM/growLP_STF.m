function htgaSTF = growLP_STF(den)

global myNTFden

myNTFden = den;

htgaOutput = HTGA(@LPSTF_obj,@LPSTF_init);

gamma = htgaOutput.F_s.best_chr;
htgaSTF.num = gamma2poly(gamma);
htgaSTF.num = htgaSTF.num/sum(htgaSTF.num)*sum(den);
    
htgaSTF.den = den;

