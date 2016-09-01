
%%
function out = cc_gamma(time,params)

n = params.exponent;
tau1 = params.tau1;
amp1 = params.amp1;
tau2 = params.tau2;
amp2 = params.amp2;

time1 = time-params.timelag1;
out1 = ((time1/tau1).^(n-1).*exp(-time1/tau1))./(tau1*factorial(n-1));
out1(time1<0) = 0;
out1 = (out1-min(out1))./ (max(out1)-min(out1));
out1 = amp1*out1;

time2 = time-params.timelag2;
out2 = ((time2/tau2).^(n-1).*exp(-time2/tau2))./(tau2*factorial(n-1));
out2(time2<0) = 0;
out2 = (out2-min(out2))./(max(out2)-min(out2));
out2 = amp2*out2;

out = out1+out2;

out = out/max(out); % scale so that peak = 1