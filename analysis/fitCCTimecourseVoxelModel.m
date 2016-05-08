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

fixedParams.disp = 0;
%% Parse and build long form stim matrix

% stimvol - basecon - newcon - basecoh - newcoh - timing
designs = {};

for di = 1:length(timeseries)
    design = zeros(10000,6);
    count = 1;
    
    cts = timeseries{di};
    csv = stimvol{di};
    crt = runtrans{di};
    cbcon = basecon{di};
    cbcoh = basecoh{di};
    cstimname = stimnames{di};
    ctim = timing{di};
    
    [con,coh,time] = parseNames(cstimname,'contrast=','coherence=','timing=',' and ');

    for ci = 1:length(csv)
        sv = csv{ci};        
        
        for si = 1:length(sv)
            design(count,:) = [sv(si) cbcon con(ci) cbcoh coh(ci) ctim(ci)];
            count = count+1;
        end
    end
    design = design(1:count,:);
    
    % internal validation
    remove_idxs = logical(logical((design(:,2)-design(:,3))==0).*logical((design(:,4)-design(:,5))==0));
    if any(remove_idxs)
        warning('You included conditions where delta con/coh == 0, removing these so the model doesn''t try to deal with them improperly');
        design = design(~remove_idxs,:);
    end
    
    designs{di} = design;
end

%% Parameter initialization
% one value = fixed
% three values = optimization
% >three values = loop across options

% Contrast Function Parameters
initparams.n = 1; 
initparams.Rmax = [1 -inf inf];
initparams.c50 = [0.1 0 1];

% Coherence Function Parameters
initparams.slope = [0.1 -inf inf];

% Constant Effect
initparams.beta = 0;

% Gamma function
initparams.amp1 = 1;
initparams.tau1 = [0.75 -inf inf];
initparams.timelag1 = [1 -inf inf];
initparams.amp2 = [-0.2 -inf 0];
initparams.tau2 = 1.2;
initparams.timelag2 = 2;
initparams.exponent = [7];
fixedParams.diff = 1;

initparams.offset = [0 -inf inf];

% Dropoff of effect
% initparams.lambda = [0.04 -inf inf];
initparams.lambda = 0;

%% Fit Model
fit = optimFitModel(initparams,timeseries,designs,runtrans);

%% Do the deconvolution analysis on the model and the fit
fit.orig.timeseries = timeseries;
fit.orig.designs = designs;
fit.orig.runtrans = runtrans;

function fit = optimFitModel(initparams,timeseries,designs,runtrans)

global fixedParams

[initparams,minparams,maxparams] = initParams(initparams);
maxiter = inf;
displsqnonlin = 'off';

optimParams = optimset('Algorithm','levenberg-marquardt','MaxIter',maxiter,'Display',displsqnonlin);

%% Do the fit

if fixedParams.disp
    f = figure;
else
    f = -inf;
end

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
        [curparams, ~, res, ~, ~, ~, curjacob] = lsqnonlin(@fitModel,initparams,minparams,maxparams,optimParams,timeseries,designs,runtrans,f);
        if res < bestres
            bestparams = curparams;
            bestres = res;
            jacobian = curjacob;
        end
    end
else
    [bestparams, ~, res, ~, ~, ~, jacobian] = lsqnonlin(@fitModel,initparams,minparams,maxparams,optimParams,timeseries,designs,runtrans,f);
end

%% Compute CIs
% jacobian = jacobian'*jacobian;
% reducedChiSquared = (res*res')/(length(hrf)-length(initparams));
% covar = sqrt(reducedChiSquared * inv(jacobian));
% 
% fit.ci_covar = nlparci(params,res,'covar',covar);
% % fit.ci_jacobian = nlparci(bestparams,res,'jacobian',jacobian);

%% Get Best Fit
[err, fit] = fitModel(bestparams,timeseries,designs,runtrans,f);
% fit.err = reshape(err,size(hrf));
fit.params = getParams(bestparams);

%% Compute full
fit.full = struct;
fcoh = 0:.1:1;
fcon = 0:.1:1;

fit.full.fcon = fcon;
fit.full.fconr = conModel(fcon,fit.params);
fit.full.fcoh = fcoh;
fit.full.fcohr = cohModel(fcoh,fit.params);

function [err, fit] = fitModel(params,timeseries,designs,runtrans,f)

global fixedParams

params = getParams(params);

% by default lets model in 0.5 s increments

t = 0:0.5:30;
impulse = gamma(t,params);

errs = {};

% For each timeseries we need to build a model timeseries. That model needs
% to not run over the transitions between runs, so we have to do some
% selection as we're building it up. This is a linear GLM using the impulse
% function above. The size of the effect at each trial is:
%
% offset + coneff + coheff
%
% and this is placed at every X time points where X is defined by the
% timing value
%
% that model timeseries is convolved with the impulse function computed
% above
%
% err is computed as timeseries - model timeseries and then simply
% concatenated across the runs
%
% stimvol - basecon - newcon - basecoh - newcoh - timing

err = zeros(length(timeseries),length(timeseries{1}));

out = cell(size(timeseries));
for ti = 1:length(timeseries)
    ts = timeseries{ti};
    ts = (ts-1)*100;
    model = zeros(size(ts));
    design = designs{ti};
    transitions = runtrans{ti};
    
    for run = 1:size(transitions,1)
        cdesign = fil(design,1,'>=',transitions(run,1));
        cdesign = fil(cdesign,1,'<=',transitions(run,2));
        
        for si = 1:size(cdesign,1)
            % okay, for each stimvol, place its effect
            sv = cdesign(si,1);
            coneff = conModel(cdesign(si,3)-cdesign(si,2),params);
            coheff = cohModel(cdesign(si,5)-cdesign(si,4),params);
            effect = coneff+coheff+params.offset;
            
            idxs = sv:min(transitions(run,2),sv+cdesign(si,6));
            model(idxs) = model(idxs)+effect;
        end
    end
    modeltimeseries = conv(impulse,model);
    modeltimeseries = modeltimeseries(1:length(ts));
    
    out{ti} = modeltimeseries;
    err(ti,:) = ts - modeltimeseries;
end

err = err(:);

fit.out = out;
fit.impulse = impulse;
% 
allts = ([timeseries{:}]-1)*100;
ssres = sum(err.^2);
sstot = sum((allts(:)-mean(allts(:))).^2);
fit.r2 = 1 - ssres/sstot;

if fixedParams.disp
    figure(f)
    clf(f)
    hold on
    plot(ts(1:1000),'r');
    plot(modeltimeseries(1:1000),'b');
    title(sprintf('R^2: %0.2f Rmax: %2.2f c50: %2.2f slope: %2.2f offset: %2.2f',fit.r2,params.Rmax,params.c50,params.slope,params.offset));
end

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