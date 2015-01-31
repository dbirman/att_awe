%% scalePDF

function [scaled, X2] = scalePdf(unscaled,X,factor)
if factor > 1
    error('Scale not designed for scale > 1');
end

X2 = min(X):factor*(X(2)-X(1)):max(X);
pad = floor((length(X2)-length(X))/2);
scaled = [zeros(1,pad) unscaled zeros(1,pad)];