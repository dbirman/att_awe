function params = fitavghrf(hrf)

time = 0.25:0.5:40.5;
% hrfparams.tau1 = [0.4 -inf inf];
% hrfparams.timelag1 = [1.5 0 3];
% hrfparams.amp2 = [-0.2 -inf 0];
% hrfparams.tau2 = [0.4 -inf inf];
% hrfparams.timelag2 = [4 0 9];

initparams = [0.4 1.5 -0.2 1 4];
minparams = [-inf 0 -inf -inf 0];
maxparams = [inf 3 0 inf 9];

[bestparams, ~, ~, ~, ~, ~, ~] = lsqnonlin(@gamRes,initparams,minparams,maxparams,[],hrf,time);

[res,params] = gamRes(bestparams,hrf,time);

function [res,hrfparams] = gamRes(params,hrf,time)


hrfparams.amp1 = 1;
hrfparams.exponent = 7;
hrfparams.tau1 = params(1);
hrfparams.timelag1 = params(2);
hrfparams.amp2 = params(3);
hrfparams.tau2 = params(4);
hrfparams.timelag2 = params(5);

fit = cc_gamma(time,hrfparams);
res = hrf-fit;