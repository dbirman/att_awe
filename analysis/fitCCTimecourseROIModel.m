function fit = fitCCTimecourseROIModel( data )
%CCROIMODEL Fit the contrast coherence model to an ROI
%
%   Dan Birman - Gardner Lab, Stanford University
%       May 10th, 2016
%
%   fit = fitCCTimecourseROIModel(data)
%
%   OUTPUT:
%   fit.r2
%   fit.params.PARAM = [-95%CI estimate +95%CI]
%   fit.full = Model output for con = 0:1 and coh 0:1
%
%   INPUT:
%   data.tSeries = tSeries from savedata_localizer
%   data.concatInfo = concatInfo from viewGet(v,'concatInfo')
%   data.design = long form design matrix [vol# basecon con basecoh coh time]
%
%   ALGORITHM:
%   We will fit a least squares model that optimizes the HRF across ROIs
%   but allows the contrast and coherence response functions to vary by
%   ROI. This seems to be the best compromise. Other ideas you tried that
%   didn't work:
%
%    - Fitting individual voxels
%    - Fitting the HRF across voxels and then fitting individual voxels
%    - Fitting to a tSeries averaged using an arbitrary R^2 cutoff
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
%   (4) The impulse response function corresponds to a 500 ms stimulus in
%   our experiments, longer stimuli cause an exponential dropoff in firing,
%   such that firing rate is modeled by:
%
%   R = Rcon(c) + Rcoh(c)
% 
%   exp(-time*lambda) * R
%

design = data.design;

if ~iscell(data.rois)
    data.rois{1} = data.rois;
end

global fixedParams
fixedParams.hrf.disp = 1;
fixedParams.vox.disp = 1;
fixedParams.diff = 1;

%% internal validation
remove_idxs = logical(logical((design(:,2)-design(:,3))==0).*logical((design(:,4)-design(:,5))==0));
if any(remove_idxs)
    warning('You included conditions where delta con/coh == 0, removing these so the model doesn''t try to deal with them improperly');
    design = design(~remove_idxs,:);
    data.design = design;
end

%% Parameter initialization
global params

% HRF Parameters
hrfparams.amp1 = 1;
hrfparams.tau1 = [0.45 -inf inf];
hrfparams.timelag1 = [.7 0 3];
hrfparams.amp2 = [-0.25 -inf 0];
hrfparams.tau2 = [1.7 -inf inf];
hrfparams.timelag2 = [0 0 6];
hrfparams.exponent = 7;
hrfparams.n = 1;

% Contrast Function Parameters
roiparams.n = 1;
roiparams.Rmax = [.1 -inf inf];
roiparams.c50 = [0.5 0 1];
% Coherence Function Parameters
roiparams.slope = [0.1 -inf inf];
% Offset
roiparams.offset = [.05 -inf inf];

params.hrfparams = hrfparams;
params.roiparams = roiparams;

params.num = length(data.ROIs);

fixedParams.sstot = 0;
for ti = 1:length(data.tSeries)
    data.tSeries{ti} = (data.tSeries{ti}-1)*100; % move into zero mean and 1% space, easier for interpretation later
    fixedParams.sstot = fixedParams.sstot + sum((data.tSeries{ti}-repmat(mean(data.tSeries{ti},2),1,size(data.tSeries{ti},2))).^2,2);
end

%% fit HRF
fit = fitModel(data);

function fit = fitModel(data)

% Fit to the mean timeseries from the top R^2 values
global params fixedParams
[hrfparams,minparams,maxparams] = initParams(params.hrfparams,'hrf');

if fixedParams.hrf.disp
    f = figure;
else
    f = -inf;
end

% Timeseries
r2mu = mean(params.r2);
r2std = std(params.r2);
r2cut = r2mu + 2*r2std;
idxs = params.sortindex(params.r2>r2cut);
while length(idxs)<50
    disp('(fitcctimecourseroimodel) Min voxels 50 not met');
    idxs = params.sortindex(params.r2>(r2cut-0.5*r2std));
end

tSeries = mean(data.tSeries(idxs,:));

optimParams = optimset('Algorithm','levenberg-marquardt','MaxIter',inf,'Display','off');
[curparams, ~, res, ~, ~, ~, curjacob] = lsqnonlin(@hrfResidual,hrfparams,minparams,maxparams,optimParams,tSeries,data.design,data.concatInfo.runTransition,'hrf',f,fixedParams);

fit.params = getParams(curparams,'hrf',fixedParams);

function res = hrfResidual(params,tSeries,design,runtrans,flag,f,fixedParams)

params = getParams(params,flag,fixedParams);

t = 0.25:0.5:49.75;
impulse = gamma(t,params);

model = zeros(size(tSeries));

for run = 1:size(runtrans,1)
    cdesign = fil(design,1,'>=',runtrans(run,1));
    cdesign = fil(cdesign,1,'<=',runtrans(run,2));

    for si = 1:size(cdesign,1)
        % okay, for each stimvol, place its effect
        sv = cdesign(si,1);
        coneff = conModel(cdesign(si,3)-cdesign(si,2),params);
        coheff = cohModel(cdesign(si,5)-cdesign(si,4),params);
        effect = coneff+coheff+params.offset;

        if cdesign(si,6)>=1
            idxs = sv:min(runtrans(run,2),sv+cdesign(si,6));
        else
            idxs = sv;
            effect = effect * cdesign(si,6);
        end
        if length(effect)>1
            model(:,idxs) = model(:,idxs)+repmat(effect,1,length(idxs));
        else
            model(:,idxs) = model(:,idxs)+repmat(effect,size(model,1),length(idxs));
        end
    end
end
cmt = conv(impulse,model);
modeltimeseries = cmt(1:size(tSeries,2));

res = tSeries-modeltimeseries;

if f>0
    figure(f)
    clf(f)
    subplot(211)
    hold on
    plot(tSeries(1:1000),'b');
    plot(modeltimeseries(1:1000),'r');
%     title(sprintf('R^2: %0.2f Rmax: %2.2f c50: %2.2f slope: %2.2f offset: %2.2f',fit.r2,params.Rmax,params.c50,params.slope,params.offset));
    subplot(212)
    plot(t,impulse);
end

ssres = sum(res.^2);

% fit.r2 = 1 - ssres/data.sstot;

%%
function out = gamma(time,params)

% global fixedParams
% 
% if fixedParams.diff
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
% else
%     n = params.exponent;
%     tau = params.tau;
%     time = time-params.timelag;
%     out = ((time/tau).^(n-1).*exp(-time/tau))./(tau*factorial(n-1));
%     out(time < 0) = 0;
%     out = (out-min(out)) ./ (max(out)-min(out));
%     out = params.amplitude*out;
% end

out = out/sum(out)/2; % normalize to sum=1 for 1% signal change / s

function out = conModel(con,params)

out = params.Rmax .* ((con.^params.n) ./ (con.^params.n + params.c50.^params.n)); 

function out = cohModel(coh,params)

out = params.slope .* coh;

function [initparams, minparams, maxparams] = initParams(params,flag)
%%
global fixedParams

fixedParams.(flag).strs = fields(params);

initparams = [];
minparams = [];
maxparams = [];
indexes = cell(1,length(fixedParams.(flag).strs));
count = 1;

fixed = zeros(1,length(fixedParams.(flag).strs));
optim = zeros(1,length(fixedParams.(flag).strs));
%%
for i = 1:length(fixedParams.(flag).strs)
    cvals = params.(fixedParams.(flag).strs{i});
    if length(cvals)==1
        fixedParams.(flag).(fixedParams.(flag).strs{i}) = cvals;
        fixed(i) = 1;
    elseif length(cvals)==3
        initparams = [initparams cvals(1)];
        minparams = [minparams cvals(2)];
        maxparams = [maxparams cvals(3)];
        indexes{i} = count;
        count = count+1;
    elseif length(cvals)==2 || length(cvals)>3
        % optimizer
        fixedParams.(flag).(fixedParams.(flag).strs{i}) = cvals;
        optim(i) = 1;
    else
        error('You initialized a parameter with the wrong initial values... unable to interpret');
    end
end
fixedParams.(flag).optim = optim;
fixedParams.(flag).fixed = fixed;
fixedParams.(flag).idx = indexes;

function p = getParams(params,flag,fixedParams)

for i = 1:length(fixedParams.(flag).strs)
    if fixedParams.(flag).fixed(i)
        p.(fixedParams.(flag).strs{i}) = fixedParams.(flag).(fixedParams.(flag).strs{i});
    else
        p.(fixedParams.(flag).strs{i}) = params(fixedParams.(flag).idx{i});
    end
end