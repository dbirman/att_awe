
function [mu, std] = buildRcurve(Y,X,bins)

binned = binData(Y,X,bins);
mu = cellfun(@nanmean,binned);
std = cellfun(@(x) sqrt(nanmean(x)*(1-nanmean(x))/length(x)),binned);