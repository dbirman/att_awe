function [params, fit] = fitCCBehavModel(adata,figs)
% FITCUMGAUSS
%
% Fit two cumulative gaussian functions (mu, std) as params using maximum
% likelihood such that the data in con/coh/resp seems most likely!
% Basically we will compute two cumulative gaussian functions and find the
% probability of every data point, then return the likelihood of that
% function.
% Mu and std describe the shapes of the gaussian functions, and alpha
% controls the proportion of probability contributed by the contrast
% function.
%             con coh
%            mu s mu s    alpha

% remove any NaN
adata = adata(~any(isnan(adata),2),:);

initparams = [0 .15 0 .15 0.5];
minparams = [-inf 0 -inf 0 0];
maxparams = [inf inf inf inf 1];

if ieNotDefined('figs')
    figs = 0;
end

% Call fmins
if figs
    f = figure;
else
    f = -1;
end
params = fmincon(@(p) fitcumgauss(p,adata,f),initparams,[],[],[],[],minparams,maxparams);

fit.params = params;
[fit.likelihood, fit.fit] = fitcumgauss(params,adata,0);

function [likelihood, fit] = fitcumgauss(params,adata,f)

con = adata(:,2);
coh = adata(:,3);
resp = adata(:,4);

likelihood = 0;

if params(2)<0, params(2)=0; end
if params(4)<0, params(4)=0; end
if params(5)<0, params(5)=0; end
if params(5)>1, params(5)=1; end

conprob = normcdf(con,params(1),params(2));
cohprob = normcdf(coh,params(3),params(4));

x = -1:.001:1;
concum = normcdf(x,params(1),params(2));
cohcum = normcdf(x,params(3),params(4));

fit.x = x;
fit.concum = concum;
fit.cohcum = cohcum;

for i = 1:length(con)
    if resp(i)
        prob = conprob(i)*params(5) + cohprob(i)*(1-params(5));
    else
        prob = (1-conprob(i))*params(5) + (1-cohprob(i))*(1-params(5));
    end
    likelihood = likelihood + log(prob);
end

likelihood = -likelihood;

if isnan(likelihood)
    keyboard
end

if f>0
    cc_rightchoice(adata,f);
end