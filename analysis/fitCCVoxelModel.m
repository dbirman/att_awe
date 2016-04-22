function fit = fitCCVoxelModel( r, con, coh )
%FITCCVOXELMOD
%
% fit = fitCCVoxelModel(amplitudes, contrasts, coherences);
%
% Dan Birman
% April 2016
%
% Fits an intercept + linear motion effect + Naka Rushton contrast effect
% to a set of data (non-linear least squares fitting). Total parameters:
%
% [Rmax c50 n slope offset]

con = con - 0.25;

% remove any no change conditions (amplitude will be 0, we don't want to
% model these)
remove_idxs = logical(logical(con==0) .* logical(coh==0));
r = r(~remove_idxs);
con = con(~remove_idxs);
coh = coh(~remove_idxs);

% % parmaeters
%              %Rmax          c50     n      slope      offset
% initParams = [max(r)        0.375     2      1       0];
% minParams =  [0             0         1      -inf     -inf];
% maxParams =  [inf           1       5      inf      inf];

             %Rmax          c50     n      slope      offset
initParams = [max(r)        0.375           1       0];
minParams =  [0             0               -inf     -inf];
maxParams =  [inf           1             inf      inf];

% optimization parameters
maxiter = inf;
optimParams = optimset('MaxIter',maxiter,'Display','off');

% now go fit
[params, ~, res, ~, ~, ~, jacobian] = lsqnonlin(@residual,initParams,minParams,maxParams,optimParams,r,con,coh);

fit = parseParams(params);

[fit.err, fit.amps] = residual(params,r,con,coh);

fit.ci = nlparci(params,res,jacobian);
fit.Rmaxint = fit.ci(1,:);
fit.c50int = fit.ci(2,:);
fit.slopeint = fit.ci(3,:);
fit.offsetint = fit.ci(4,:);

x = 0:.01:1;
fit.range = x;
fit.coneff = nakaRushton(x,fit);
fit.coheff = fit.slope*x;

ssres = sum(fit.err.^2);
sstot = sum((r-mean(r)).^2);
fit.r2 = 1- ssres/sstot;

fcoh = [0 1];
fcon = [repmat(0,1,length(fcoh)) repmat(0.25,1,length(fcoh)) repmat(0.5,1,length(fcoh)) repmat(0.75,1,length(fcoh))];

fcoh = repmat(fcoh,1,4);

fit.full.con = fcon;
fit.full.coh = fcoh;
fit.full.amps = fit.offset + fit.slope*fcoh + nakaRushton(fcon,fit);
fit.full.conamps = nakaRushton(fcon,fit);
fit.full.cohamps = fit.slope * fcoh;

function [err, fit] = residual(params,r,con,coh)

p = parseParams(params);

% intercept + linear effect of con + naka effect of contrast
fit = p.offset + p.slope * coh + nakaRushton(con,p);

err = fit - r;

%%%%%%%%%%%%%%%%%%%%%
%    nakaRushton    %
%%%%%%%%%%%%%%%%%%%%%
function response = nakaRushton(c,p)

% response = p.Rmax * ((c.^p.n) ./ ((c.^p.n) + p.c50.^p.n));

response = (p.Rmax * c) ./ (p.c50 + c);


%%%%%%%%%%%%%%%%%%%%%
%    parseParams    %
%%%%%%%%%%%%%%%%%%%%%
function p = parseParams(params)

% if m.fixedN
%   p.Rmax = params(1);
%   p.c50 = params(2);
%   p.n = m.n;
%   p.offset = params(3);
% else
  p.Rmax = params(1);
  p.c50 = params(2);
%   p.n = params(3);
  p.slope = params(3);
  p.offset = params(4);
% end  


