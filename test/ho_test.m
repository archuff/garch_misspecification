function results = ho_test(e,h,ord,param)

p = ord(2);
q = ord(1);
e2 = e.^2;
T = length(e2);

dB = [e e.^3];
sres = e2./h - 1;
dgarch = [ones(T,1) e2 h];

stat = ho_test_core(T,h,param,dB,dgarch,sres);

pval = 1-chi2cdf(stat,2);
results.pval = pval;
results.stat = stat;

end
