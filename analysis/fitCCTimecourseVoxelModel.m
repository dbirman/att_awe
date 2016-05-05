function fit = fitCCTimecourseVoxelModel( timeseries, stimvol, runtrans, basecon, basecoh, stimnames, timing )
%CCVOXELMODEL Fit the contrast coherence model to a voxel or average of
%voxels
%
%   Dan Birman - Gardner Lab, Stanford University
%       April 20,2016
%
%   fit = ccVoxelModel(hrfData, con, coh)
%
%   OUTPUT:
%   fit.r2
%   fit.params.PARAM = [-95%CI estimate +95%CI]
%   fit.full = Model output for con = 0:1 and coh 0:1
%
%   INPUT:
%   timecourse = each cell is a timecourse that will be fit
%   stimvol = each cell is a set of stimvols that will be used
%   bases = each cell is the base con/coh for those conditions
%   con/coh = each cell is the conditions corresponding to the stimvols
%   timing = each cell contains for each condition the stimulus length (in
%   250 ms volumes)
%
%   Based on Rees & Koch 2000 and a bunch of other stuff...
%
%   Assumptions:
%   (1) Contrast influences firing rates in a non-linear manner, according
%   to the function:
%
%   Rcon(c) = Rmax * c^n / (c^n + c50^n)
%
%   (2) Coherence influences firing rates in a linear manner, according to
%   the function:
%
%   Rcoh(c) = beta * c
%
%   (3) The amplitude of the human fMRI BOLD response is a linear function 
%   of the rate of spiking of neurons within each voxel. Therefore the
%   impulse response function due to a particular contrast and coherence:
%
%   HRF(con,coh) = conv(Rcon(c) + Rcoh(c) + offset,Canonical)
%
%   Where Canonical is a gamma function with three parameters: exponent,
%   tau, and timelag.
%
%   (4) The impulse response function corresponds to a 250 ms stimulus in
%   our experiments, longer stimuli cause an exponential dropoff in firing,
%   such that firing rate is modeled by:
%
%   R = Rcon(c) + Rcoh(c)
%
%   exp(-time*lambda) * R

% the fixedParams will help keep track of parameters that are fixed
global fixedParams

%% Deal with multiple timeseries constraint satisfaction



%% Validation
if max(con)>1 || min(con)<0 || max(coh)>1 || min(coh)<0 || size(hrf,1) ~= size(con,1) || size(hrf,1) ~= size(coh,1) || size(hrf,1) ~= size(timing,1) || size(hrf,1) ~= size(basecon,1) || size(hrf,1) ~= size(basecoh,1)
    error('Data improperly formatted');
end

remove_idxs = logical(logical((basecon-con)==0) .* logical((basecoh-coh)==0));
if any(remove_idxs)
    warning('You included conditions where con=0 and coh=0 (i.e. no change occurred), the model can''t account for these so they will be removed');  
    hrf = hrf(~remove_idxs,:);
    con = con(~remove_idxs);
    coh = coh(~remove_idxs);
    basecon = basecon(~remove_idxs);
    basecoh = basecoh(~remove_idxs);
    timing = timing(~remove_idxs);
end

%% Parameter initialization
% one value = fixed
% three values = optimization
% >three values = loop across options

% Contrast Function Parameters
if ~any((con-basecon)>0)
    % fix contrast
    warning('No contrast changes: Fixing to 0 effect');
    initparams.n = 1; initparams.Rmax = 0; initparams.c50 = 0.5;
else
    initparams.n = 1; % one value for fixed
    initparams.Rmax = [1 0 inf]; % three values for optimization
    initparams.c50 = [0.1 -inf inf];
end

% Coherence Function Parameters
initparams.slope = [0.1 -inf inf];

% Constant Effect
initparams.beta = 0;

% Gamma function
initparams.amp1 = 1;
initparams.tau1 = [0.75 0.5 1];
initparams.timelag1 = 0.75;
initparams.amp2 = 0;
initparams.tau2 = 1.2;
initparams.timelag2 = 2;
initparams.exponent = 6;
fixedParams.diff = 1;

% initparams.amplitude = 0;
% initparams.tau = 0;
% initparams.timelag = 0;

initparams.offset = [0 0 inf];

% Dropoff of effect
% initparams.lambda = [0.04 -inf inf];
initparams.lambda = 0;

%% Parameter Fixing
if ~any((coh-basecoh)>0)
    % fix coherence
    warning('No coherence changes: Fixing to 0 effect');
    initparams.slope = 0;
end

%% Fit Model
fit = optimFitModel(initparams,hrf,basecon,basecoh,con,coh,timing,time,runtrans);

fit.orig.hrf = hrf;
fit.orig.con = con;
fit.orig.coh = coh;
fit.orig.timing = timing;
fit.orig.basecon = basecon;
fit.orig.basecoh = basecoh;
fit.orig.time = time;

function fit = optimFitModel(initparams,hrf,basecon,basecoh,con,coh,timing,time)

global fixedParams

[initparams,minparams,maxparams] = initParams(initparams);
maxiter = inf;
displsqnonlin = 'off';

optimParams = optimset('Algorithm','levenberg-marquardt','MaxIter',maxiter,'Display',displsqnonlin);

%% Do the fit

if any(fixedParams.optim)
    % Deal with the parameters that need to be optimized
    bestparams = struct;
    bestres = inf;
    jacobian = [];
    if sum(fixedParams.optim)>1, error('Can''t optimize over multiple parameters'); end
    idx = find(fixedParams.optim,1);
    cpvals = fixedParams.(fixedParams.strs{idx});
    for i = 1:length(cpvals)
        fixedParams.(fixedParams.strs{idx}) = cpvals(i);
        fixedParams.fixed(idx) = 1;
        [curparams, ~, res, ~, ~, ~, curjacob] = lsqnonlin(@fitModel,initparams,minparams,maxparams,optimParams,hrf,basecon,basecoh,con,coh,timing,time);
        if res < bestres
            bestparams = curparams;
            bestres = res;
            jacobian = curjacob;
        end
    end
else
    [bestparams, ~, res, ~, ~, ~, jacobian] = lsqnonlin(@fitModel,initparams,minparams,maxparams,optimParams,hrf,basecon,basecoh,con,coh,timing,time);
end

%% Compute CIs
% jacobian = jacobian'*jacobian;
% reducedChiSquared = (res*res')/(length(hrf)-length(initparams));
% covar = sqrt(reducedChiSquared * inv(jacobian));
% 
% fit.ci_covar = nlparci(params,res,'covar',covar);
% % fit.ci_jacobian = nlparci(bestparams,res,'jacobian',jacobian);

%% Get Best Fit
[err, fit] = fitModel(bestparams,hrf,basecon,basecoh,con,coh,timing,time);
fit.err = reshape(err,size(hrf));
fit.params = getParams(bestparams);

%% Compute full
fit.full = struct;
fcoh = 0:.1:1;
fcon = 0:.1:1;

fit.full.fcon = fcon;
fit.full.fconr = conModel(fcon,fit.params);
fit.full.fcoh = fcoh;
fit.full.fcohr = cohModel(fcoh,fit.params);

function [err, fit] = fitModel(params,hrf,basecon,basecoh,con,coh,convtiming,t)
params = getParams(params);

impulse = gamma(t,params);

out = zeros(size(hrf));
err = zeros(size(hrf));
for hi = 1:size(hrf,1)
    % get the size of the 'neural firing rate change'
    coneff = conModel(con(hi),params)-conModel(basecon(hi),params);
    coheff = cohModel(coh(hi),params)-cohModel(basecoh(hi),params);
    % total effect
    effect = coneff + coheff + params.offset;
    % timing
    timing = zeros(size(impulse));
%     timing(1:convtiming(hi)) = exp(-(1:convtiming(hi))*params.lambda);
    timing(1:convtiming(hi)) = ones(1,convtiming(hi))*params.lambda + effect;
    timing(1) = timing(1) + params.beta;
    % convolve
    out_ = conv(impulse,timing);
    out(hi,:) = out_(1:size(hrf,2));
    err(hi,:) = hrf(hi,:)-out(hi,:);
end

% err = hrf-out;
err = err(:);

fit.out = out;
fit.impulse = impulse;

ssres = sum(err.^2);
sstot = sum((hrf(:)-mean(hrf(:))).^2);
fit.r2 = 1 - ssres/sstot;

%%
function out = gamma(time,params)

global fixedParams

if fixedParams.diff
    n = params.exponent;
    tau1 = params.tau1;
    amp1 = params.amp1;
    tau2 = params.tau2;
    amp2 = params.amp2;
    
    time1 = time-params.timelag1;
    out1 = ((time1/tau1).^(n-1).*exp(-time1/tau1))./(tau2*factorial(n-1));
    out1(time1<0) = 0;
    out1 = (out1-min(out1))./ (max(out1)-min(out1));
    out1 = amp1*out1;
    
    time2 = time-params.timelag2;
    out2 = ((time2/tau2).^(n-1).*exp(-time2/tau2))./(tau2*factorial(n-1));
    out2(time2<0) = 0;
    out2 = (out2-min(out2))./(max(out2)-min(out2));
    out2 = amp2*out2;
    
    out = out1+out2;
else
    n = params.exponent;
    tau = params.tau;
    time = time-params.timelag;
    out = ((time/tau).^(n-1).*exp(-time/tau))./(tau*factorial(n-1));
    out(time < 0) = 0;
    out = (out-min(out)) ./ (max(out)-min(out));
    out = params.amplitude*out;
end

function out = conModel(con,params)

out = params.Rmax .* ((con.^params.n) ./ (con.^params.n + params.c50.^params.n)); 

function out = cohModel(coh,params)

out = params.slope .* coh;

function [initparams, minparams, maxparams] = initParams(params)

global fixedParams

fixedParams.strs = fields(params);
fixedParams.num = length(fixedParams.strs);

initparams = [];
minparams = [];
maxparams = [];
indexes = zeros(1,fixedParams.num);
count = 1;

fixed = zeros(1,fixedParams.num);
optim = zeros(1,fixedParams.num);

for i = 1:fixedParams.num
    cvals = params.(fixedParams.strs{i});
    
    if length(cvals)==1
        fixedParams.(fixedParams.strs{i}) = cvals;
        fixed(i) = 1;
    elseif length(cvals)==3
        initparams = [initparams cvals(1)];
        minparams = [minparams cvals(2)];
        maxparams = [maxparams cvals(3)];
        indexes(i) = count;
        count = count+1;
    elseif length(cvals)==2 || length(cvals)>3
        % optimizer
        fixedParams.(fixedParams.strs{i}) = cvals;
        optim(i) = 1;
    else
        error('You initialized a parameter with the wrong initial values... unable to interpret');
    end
end
fixedParams.optim = optim;
fixedParams.fixed = fixed;
fixedParams.idx = indexes;

function p = getParams(params)

global fixedParams

for i = 1:fixedParams.num
    if fixedParams.fixed(i)
        p.(fixedParams.strs{i}) = fixedParams.(fixedParams.strs{i});
    else
        p.(fixedParams.strs{i}) = params(fixedParams.idx(i));
    end
end