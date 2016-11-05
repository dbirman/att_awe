
function [mu, std] = buildRcurve(Y,X,bins)

binned = binData(Y,X,bins);
for bi = 1:length(binned)
    if length(binned{bi})<20
        binned{bi} = [];
    end
end
% cis = cellfun(@(x) bootci(10000,@nanmean,x),binned);
% mu = nanmean(binned);
mu = cellfun(@nanmean,binned);
% std = cis(2,:)-mu;
std = cellfun(@(x) sqrt(nanmean(x)*(1-nanmean(x))/length(x)),binned);
