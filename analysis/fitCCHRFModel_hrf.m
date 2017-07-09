function fit = fitCCHRFModel_hrf( data, mode, dataopt)
%CCROIMODEL Fit the contrast coherence model to an ROI
%
%   Dan Birman - Gardner Lab, Stanford University
%       May 10th, 2016
%
%   fit = fitCCHRFModel(data)
%
%   OUTPUT:

%% Setup

global fixedParams
fixedParams = struct;
fixedParams.ROIs = data.ROIs;

%% Set data
data.cc.cresp = data.cc.(dataopt);
data.time.cresp = data.time.(dataopt);

%% parse mode:
hrfparams = struct;
roiparams = struct;

if strfind(mode,'spkdec')
    % Spike rate decay
    hrfparams.spkexp = -0.623;%[-0.5 -inf 0];
    hrfparams.hrfexp = 0;
    fixedParams.spkdec = 1;
elseif strfind(mode,'spkhrfdec')
    hrfparams.spkexp = -0.623;
    hrfparams.hrfexp = -0.623;
else
    % HRF response decay
    hrfparams.spkexp = 0;
    hrfparams.hrfexp = [-0.623 -inf inf]; % adaptation exponent (for time)
end

for i = 1:(8*(20+size(data.time.cresp,2)))
    hrfparams.(sprintf('betas%i',i)) = [1 -inf inf];
end

hrfparams.amp1 = 1;
hrfparams.tau1 = [0.4 -inf inf];
hrfparams.timelag1 = 1.5;
hrfparams.amp2 = [-0.3 -inf 0];
hrfparams.tau2 = [0.6 -inf inf];
hrfparams.timelag2 = 4;
hrfparams.exponent = 7;

fixedParams.regularize=0;
if strfind(mode,'doreg')
    fixedParams.regularize = 1;
end

%% Parameter initialization
global params

params.hrfparams = hrfparams;
params.roiparams = roiparams;
% 
adat = [data.cc.cresp(:); data.time.cresp(:)];
fixedParams.sstot = sum((adat-mean(adat)).^2);

%% fit HRF
fit = fitModel(data);

fit.mode = mode;
if isfield(data,'basecon')
    fit.basecon = data.basecon;
    fit.basecoh = data.basecoh;
end

function fit = fitModel(data)
%%
% Fit to the mean timeseries from the top R^2 values
global fixedParams
[initparams,minparams,maxparams] = initParams;

f = figure;

optimParams = optimset('Algorithm','trust-region-reflective','MaxIter',inf,'Display','off');
[bestparams, ~, ~, ~, ~, ~, ~] = lsqnonlin(@hrfResidual,initparams,minparams,maxparams,optimParams,data,f,fixedParams);

n = length(data.cc.cresp(:))+length(data.time.cresp(:));

[~,fit] = hrfResidual(bestparams,data,-1,fixedParams);
% 
% fit.r = corrcoef([[fit.cc.cresp(:); fit.time.cresp(:)] [fit.cc.model(:); fit.time.model(:)]]);
% fit.r2 = (fit.r(1,2))^2;

% fit.BIC = n*log(fit.SSE/n) + length(bestparams)*log(n);
% fit.AIC = n*log(fit.SSE/n) + length(bestparams)*2;
% fit.like = log(fit.SSE/n);
fit.params = getParams(bestparams,fixedParams);
% fit.params = getROIParams(fit.params,data.ROIs{1});
fit.ROIs = fixedParams.ROIs;

function [res, fit] = hrfResidual(params,data,f,fixedParams)
%%
fit = struct;

% stimvol basecon lcon rcon basecoh lcoh rcoh timing task
params = getParams(params,fixedParams);

% use non-canonical    
t = 0.25:0.5:40.5;
impulse = cc_gamma(t,params);

fit.impulse = impulse;

data.utimes = [0.5 1 2 4 5 8];
data.reps = [1 1 2 4 5 8];

for ui = 1:length(data.utimes)
    events = repmat(data.utimes(ui)^params.hrfexp,1,data.reps(ui));
    if ui==1, events = events*0.5; end
    canon = conv(events,fit.impulse);
    data.canonical(ui,:) = canon(1:length(fit.impulse));
end

% fit.cc = data.cc;
% fit.time = data.time;
% fit.cc.model = zeros(size(fit.cc.cresp));
% fit.time.model = zeros(size(fit.time.cresp));

res = zeros(1,size(data.canonical,2)*8*40);

cc_model = zeros(size(data.cc.cresp));
cc_res = zeros(size(data.cc.cresp));

for ri = 1:8
    for gi = 1:size(cc_model,2)
        beta = params.(sprintf('betas%i',(ri-1)*20+gi));
        sidx = ((ri-1)*20+gi-1)*81+1;
        cc_model(ri,gi,:) = data.canonical(data.cc.time(gi)==data.utimes,:)*beta;
        cc_res(ri,gi,:) = cc_model(ri,gi,:) - data.cc.cresp(ri,gi,:);
    end
end

time_model = zeros(size(data.time.cresp));
time_res = zeros(size(data.time.cresp));

for ri = 1:8
    for gi = 1:size(time_model,2)
        beta = params.(sprintf('betas%i',(ri-1)*size(time_model,2)+gi));
        sidx = ((ri-1)*(size(time_model,2))+gi-1)*81+1;
        time_model(ri,gi,:) = data.canonical(data.time.time(gi)==data.utimes,:)*beta;
        time_res(ri,gi,:) = time_model(ri,gi,:) - data.time.cresp(ri,gi,:);
    end
end


fit.cc.model = cc_model;
fit.time.model = time_model;

res = [cc_res(:)' time_res(:)'];
 
fit.y = [cc_model(:);time_model(:)];
fit.y_ = [data.cc.cresp(:); data.time.cresp(:)];
fit.r2 = myr2(fit.y,fit.y_);

if f>0
    figure(f)
    subplot(2,2,1:2);
    plot(res);
%     subplot(3,2,5);
%     hold on
%     plot(squeeze(fit.cc.cresp)');
%     plot(squeeze(cc_model)');
    subplot(2,2,3:4);
    plot(fit.impulse);
    
    title(sprintf('R^2: %01.3f%%',fit.r2*100));
end
% 
% fit.likelihood = ssres;

function [initparams, minparams, maxparams] = initParams()
%%
global fixedParams params

%% Deal with HRF params
fixedParams.strs = fields(params.hrfparams)';

initparams = [];
minparams = [];
maxparams = [];
indexes = {};
count = 1;

fixed = [];
optim = zeros(size(indexes));

count = 1;
for i = 1:length(fixedParams.strs)
    cvals = params.hrfparams.(fixedParams.strs{i});
    if length(cvals)==1
        fixedParams.(fixedParams.strs{i}) = cvals;
        fixed(end+1) = 1;
        indexes{end+1} = [];
    elseif length(cvals)==3
        initparams(end+1) = cvals(1);
        minparams(end+1) = cvals(2);
        maxparams(end+1) = cvals(3);
        fixed(end+1) = 0;
        indexes{end+1} = count; count = count+1;
    elseif length(cvals)==2 || length(cvals)>3
        warning('Failure');
        keyboard
        % optimizer
        fixedParams.(fixedParams.strs{i}) = cvals;
        optim(i) = 1;
    else
        error('You initialized a parameter with the wrong initial values... unable to interpret');
    end
end

%% Deal with ROI params
if iscell(params.roiparams)
    rStrs = {};
    disp('ROI parameters were passed in as a cell: interpreting as fixed values');
    for ri = 1:length(fixedParams.ROIs)
        % get the fields
        crfields = fields(params.roiparams{ri});
        % copy each field with this ROIs prefix
        for ci = 1:length(crfields)
            rStrs{end+1} = sprintf('%s%s',fixedParams.ROIs{ri},crfields{ci});
            fixedParams.(rStrs{end}) = params.roiparams{ri}.(crfields{ci});
            fixed(end+1) = 1;
            indexes{end+1} = [];
        end
    end
else
    rfields = fields(params.roiparams);

    rStrs = {};
    for ni = 1:length(rfields)
        cvals = params.roiparams.(rfields{ni});
        % check if this name includes an ROI already, if so, just copy it in
        % directly
        cfield = rfields{ni};
        if any(cellfun(@(x) ~isempty(strfind(cfield,x)),fixedParams.ROIs))
            % just copy directly
            rStrs{end+1} = cfield;
            fixedParams.(rStrs{end}) = cvals;
            fixed(end+1) = 1;
            indexes{end+1} = [];
        else
            % replicate for every ROI
            for ri = 1:length(fixedParams.ROIs)
                rStrs{end+1} = sprintf('%s%s',fixedParams.ROIs{ri},rfields{ni});
                if length(cvals)==1
                    fixedParams.(rStrs{end}) = cvals;
                    fixed(end+1) = 1;
                    indexes{end+1} = [];
                elseif length(cvals)==3
                    initparams(end+1) = cvals(1);
                    minparams(end+1) = cvals(2);
                    maxparams(end+1) = cvals(3);
                    indexes{end+1} = count;
                    count = count+1;
                    fixed(end+1) = 0;
                elseif length(cvals)==2 || length(cvals)>3
                    error('You are not allowed to use the optimizer for ROI specific parameters...');
                else
                    error('You initialized a parameter with the wrong initial values... unable to interpret');
                end
            end
        end
    end 
end

%% Save optim/fixed/indexes
fixedParams.strs = [fixedParams.strs rStrs];
fixedParams.optim = optim;
fixedParams.fixed = fixed;
fixedParams.idx = indexes;

function p = getROIParams(params,ROI)
% just grab anything that starts with ROI
p = struct;
flds = fields(params);
for fi = 1:length(flds)
    if strfind(flds{fi},ROI)
        p.(strrep(flds{fi},ROI,'')) = params.(flds{fi});
    end
end

function p = getParams(params,fixedParams)

for i = 1:length(fixedParams.strs)
    if fixedParams.fixed(i)
        p.(fixedParams.strs{i}) = fixedParams.(fixedParams.strs{i});
    else
        p.(fixedParams.strs{i}) = params(fixedParams.idx{i});
    end
end