function fit = fitCCHRFModel( data , mode)
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

%% If multiple ROIs, fit each individually
if length(data.ROIs)>1
    for ri = 1:length(data.ROIs)
        disp(sprintf('(fitCCHRFModel) Running ROI: %s',data.ROIs{ri}));
        
        ndata = data;
        ndata.ROIs = ndata.ROIs(ri);
        ndata.cc.resp = ndata.cc.resp(ri,:,:);
        ndata.time.resp = ndata.time.resp(ri,:,:);
        
        fit.roifit{ri} = fitCCHRFModel(ndata,mode);
        fit.r2(ri) = fit.roifit{ri}.r2;
        fit.BIC(ri) = fit.roifit{ri}.BIC;
    end
   
    return
end
%% parse mode:
hrfparams = struct;
roiparams = struct;
fixedParams.fitroi = 0;

hrfparams.exponent = [-0.7]; % adaptation exponent (for time)

fixedParams.numparams = 0;
if strfind(mode,'fitroi')
    % get hrf params
    % Offset
    if isempty(strfind(mode,'nooffset'))
        roiparams.offset = [0 -inf inf];
        fixedParams.numparams = fixedParams.numparams+1;
    else
        roiparams.offset=0;
    end
    % Contrast Function Parameters
    roiparams.conalpha = [1 -inf inf];
    roiparams.conkappa = [1 -inf inf];
    roiparams.conmodel = 3;
    % Coherence Function Parameters
    roiparams.cohalpha = [1 -inf inf];
    roiparams.cohkappa = [1 -inf inf];
    roiparams.cohmodel = 3;
    fixedParams.numparams = fixedParams.numparams+4;
    % run type
    fixedParams.fitroi = 1;
end

fixedParams.regularize=0;
if strfind(mode,'noreg')
    fixedParams.regularize = 1;
end

%% Parameter initialization
global params

params.hrfparams = hrfparams;
params.roiparams = roiparams;
% 
adat = [data.cc.resp(:); data.time.resp(:)];
fixedParams.sstot = sum((adat-mean(adat)).^2);

%% Change base contrast if >0
if data.basecon>0
    data.realbasecon = data.basecon;
    data.cc.con = data.cc.con-data.basecon;
    data.time.con = data.time.con-data.basecon;
    data.basecon = 0;
end
%% fit HRF
fit = fitModel(data);

fit.mode = mode;
if isfield(data,'basecon')
    fit.basecon = data.basecon;
    fit.basecoh = data.basecoh;
end

function fit = fitModel(data)

% Fit to the mean timeseries from the top R^2 values
global fixedParams
[initparams,minparams,maxparams] = initParams;

data.utimes = [0.5 1 2 4 5 8];
data.reps = [1 1 2 4 5 8];

for ui = 1:length(data.utimes)
    events = repmat(data.utimes(ui)^fixedParams.exponent,1,data.reps(ui));
    if ui==1, events = events*0.5; end
    canon = conv(events,data.hrf);
    data.canonical(ui,:) = canon(1:length(data.hrf));
end

% Transform data to fitted version
out_cc = zeros(1,size(data.cc.resp,2));
for i=1:size(data.cc.resp,2)
    curhrf = data.canonical(data.utimes==data.cc.time(i),:);
    out_cc(i) = curhrf'\squeeze(data.cc.resp(1,i,:));
end
data.cc.resp_ = out_cc;

out_time = zeros(1,size(data.time.resp,2));
for i=1:size(data.time.time,2)
    curhrf = data.canonical(data.utimes==data.time.time(i),:);
    out_time(i) = curhrf'\squeeze(data.time.resp(1,i,:));
end
data.time.resp_ = out_time;

f = figure;

optimParams = optimset('Algorithm','trust-region-reflective','MaxIter',inf,'Display','off');
[bestparams, ~, ~, ~, ~, ~, ~] = lsqnonlin(@hrfResidual,initparams,minparams,maxparams,optimParams,data,f,fixedParams);

n = length(data.cc.resp(:))+length(data.time.resp(:));

[~,fit] = hrfResidual(bestparams,data,-1,fixedParams);
fit.SSE = sum(fit.rres.^2);
fit.sstot = fixedParams.sstot;
fit.r2 = 1 - (fit.SSE/fit.sstot);
fit.BIC = n*log(fit.SSE/n) + fixedParams.numparams*log(n);
fit.params = getParams(bestparams,fixedParams);
fit.params = getROIParams(fit.params,data.ROIs{1});
fit.ROIs = fixedParams.ROIs;

function [res, fit] = hrfResidual(params,data,f,fixedParams)

fit = struct;

% stimvol basecon lcon rcon basecoh lcoh rcoh timing task
params = getParams(params,fixedParams);

fit.cc = data.cc;
fit.time = data.time;
fit.cc.model = zeros(size(fit.cc.resp));
fit.time.model = zeros(size(fit.time.resp));

res = zeros(1,length(data.cc.resp_)+length(data.time.resp_));

roiparams = getROIParams(params,fixedParams.ROIs{1});

baseConResp = conModel(data.basecon,roiparams,0,1);
baseCohResp = cohModel(data.basecoh,roiparams,0,1);

cc_model = zeros(size(data.cc.resp));

for i = 1:length(data.cc.resp_)
    ccon = data.cc.con(i);
    ccoh = data.cc.coh(i);
    
    conEff = conModel(ccon,roiparams,0,1)-baseConResp;
    cohEff = cohModel(ccoh,roiparams,0,1)-baseCohResp;
    
    effect = conEff+cohEff+roiparams.offset;
    
    cc_model(1,i,:) = data.canonical(data.cc.time(i)==data.utimes,:)*effect;
    
    res(i) = effect-data.cc.resp_(i);
end

time_model = zeros(size(data.time.resp));

for i = 1:length(data.time.resp_)
    ccon = data.time.con(i);
    ccoh = data.time.coh(i);
    
    conEff = conModel(ccon,roiparams,0,1)-baseConResp;
    cohEff = cohModel(ccoh,roiparams,0,1)-baseCohResp;
    
    effect = conEff+cohEff+roiparams.offset;
    
    time_model(1,i,:) = data.canonical(data.time.time(i)==data.utimes,:)*effect;
    
    res(length(data.cc.resp_)+i) = effect-data.time.resp_(i);
end

if fixedParams.regularize
    res = [res 0.1*conModel(0:.1:1,roiparams,0,1) 0.1*cohModel(0:.1:1,roiparams,0,1)];
end
    
if f>0
    figure(f)
    subplot(2,2,1:2);
    plot(res);
    subplot(2,2,3);
    x = 0:.01:1;
    plot(x,conModel(x,roiparams,0,1));
    subplot(2,2,4);
    plot(x,cohModel(x,roiparams,0,1));
end

fit.cc.model = cc_model;
fit.time.model = time_model;
 
try
    rres = [cc_model(:); time_model(:)] - [data.cc.resp(:); data.time.resp(:)];
catch
    keyboard
end
ssres = sum(rres.^2);
fit.rres = rres;
fit.r2 = 1 - ssres/fixedParams.sstot;
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