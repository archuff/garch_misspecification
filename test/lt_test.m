function results = lt_test(e,h,param,ord)

p = ord(1);
q = ord(2);	
nblag = p + q + 1;
T = length(e);
n = 1;
e2 = e.^2;        
sres = e2./h - 1;
dof = (n+1)*q;
dB = [e e.^3];
dG = [ones(T,1) e2 h];
stat = lt_test_core(T,h,param,dB,dG,sres,nblag);

results.stat = stat;
results.pval = 1-chi2cdf(stat,dof);  
results.dof = dof;

end
