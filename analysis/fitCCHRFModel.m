function fit = fitCCHRFModel( data , mode, pfit)
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
        if strfind(mode,'fitexp')
            ndata.params = pfit.roifit{ri}.params;
        end
        if strfind(mode,'fitsigma')
            fit.roifit{ri} = fitCCHRFModel(ndata,mode,pfit);
        else
            fit.roifit{ri} = fitCCHRFModel(ndata,mode);
        end
        fit.r2(ri) = fit.roifit{ri}.r2;
        fit.BIC(ri) = fit.roifit{ri}.BIC;
        fit.AIC(ri) = fit.roifit{ri}.AIC;
        fit.like(ri) = fit.roifit{ri}.like;
    end
   
    return
end
%% parse mode:
hrfparams = struct;
roiparams = struct;
fixedParams.fitroi = 0;
fixedParams.spkdec = 0;
fixedParams.fitexp = 0;

if strfind(mode,'fitsigma')
    if strfind(mode,'nooffset')
        roiparams.offset=0;
    elseif strfind(mode,'doubleoffset')
        roiparams.conoffset = [0 -inf inf];
        roiparams.cohoffset = [0 -inf inf];
    else
        roiparams.offset = [0 -inf inf];
    end
    fixedParams.x = pfit.x;
    fixedParams.con = pfit.con;
    fixedParams.coh = pfit.coh;
    roiparams.conmodel = 4;
    roiparams.cohmodel = 4;
    roiparams.sigmacon = [0.1 0 1];
    roiparams.sigmacoh = [0.1 0 1];
else
    roiparams.sigmacon = 1;
    roiparams.sigmacoh = 1;
end

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
    hrfparams.hrfexp = -0.623; % adaptation exponent (for time)
end

fixedParams.numparams = 0;
if strfind(mode,'fitroi')
    roiparams.conRmax = [2 -inf inf];
    roiparams.conp = 0.3;%[ eps inf];
    roiparams.conq = 2;
    roiparams.conc50 = [0.5 eps 1-eps];
    roiparams.conmodel = 2;
    % get hrf params
    % Offset
    
    if strfind(mode,'nooffset')
        roiparams.offset=0;
    elseif strfind(mode,'doubleoffset')
        roiparams.conoffset = [0 -inf inf];
        roiparams.cohoffset = [0 -inf inf];
    else
        roiparams.offset = [0 -inf inf];
    end
    
    if strfind(mode,'linear')
%         roiparams.conslope = [1 -inf inf];
%         roiparams.conmodel = 1;
        roiparams.cohslope = [1 -inf inf];
        roiparams.cohmodel = 1;
    elseif strfind(mode,'naka')
        roiparams.cohRmax = [2 -inf inf];
        roiparams.cohp = 0.3;%[2 eps inf];
        roiparams.cohq = 2;
        roiparams.cohc50 = [0.5 eps 1-eps];
        roiparams.cohmodel = 2;
    elseif strfind(mode,'tanh')
        roiparams.conalpha = [1 -inf inf];
        roiparams.conkappa = [1 -inf inf];
        roiparams.conmodel = 5;
        roiparams.cohalpha = [1 -inf inf];
        roiparams.cohkappa = [1 -inf inf];
        roiparams.cohmodel = 5;
    elseif strfind(mode,'explin')
        roiparams.conalpha = [1 0 inf];
        roiparams.conkappa = [1 0 inf];
        roiparams.conmodel = 6;
        roiparams.cohalpha = [1 0 inf];
        roiparams.cohkappa = [1 0 inf];
        roiparams.cohmodel = 6;
    elseif strfind(mode,'exp')
        % Contrast Function Parameters
        roiparams.conalpha = [1 -inf inf];
        roiparams.conkappa = [1 -inf inf];
        roiparams.conmodel = 3;
        % Coherence Function Parameters
        roiparams.cohalpha = [1 -inf inf];
        roiparams.cohkappa = [1 -inf inf];
        roiparams.cohmodel = 3;
    else
        disp('No model specified');
        keyboard
    end
    % run type
    fixedParams.fitroi = 1;
end

if strfind(mode,'interaction')
    roiparams.inbeta = [0.1 -inf inf];
else
    roiparams.inbeta = 0;
end

fixedParams.regularize=0;
if strfind(mode,'doreg')
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
%%
% Fit to the mean timeseries from the top R^2 values
global fixedParams
[initparams,minparams,maxparams] = initParams;

data.utimes = [0.5 1 2 4 5 8];
data.reps = [1 1 2 4 5 8];

if fixedParams.fitexp
    data.canonical = data.hrf;
else
    for ui = 1:length(data.utimes)
        events = repmat(data.utimes(ui)^fixedParams.hrfexp,1,data.reps(ui));
        if ui==1, events = events*0.5; end
        canon = conv(events,data.hrf);
        data.canonical(ui,:) = canon(1:length(data.hrf));
    end
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
fit.BIC = n*log(fit.SSE/n) + length(bestparams)*log(n);
fit.AIC = n*log(fit.SSE/n) + length(bestparams)*2;
fit.like = log(fit.SSE/n);
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

baseConResp = conModel(data.basecon,roiparams);
baseCohResp = cohModel(data.basecoh,roiparams);

cc_model = zeros(size(data.cc.resp));

for i = 1:length(data.cc.resp_)
    ccon = data.cc.con(i);
    ccoh = data.cc.coh(i);
    
    conEff = roiparams.sigmacon*(conModel(ccon,roiparams)-baseConResp);
    cohEff = roiparams.sigmacoh*(cohModel(ccoh,roiparams)-baseCohResp);
    inEff = roiparams.inbeta*conEff*cohEff;
    
    if conEff==0 && cohEff==0 % no change! res=0
        res(i) = 0;
    else
        if isfield(roiparams,'offset')
            if (conEff>0) || (cohEff>0)
                effect = conEff+cohEff+roiparams.offset;
            else
                effect = conEff+cohEff;
            end
        else
            effect = conEff+cohEff;
            if conEff>0
                effect = effect+ roiparams.conoffset;
            end
            if cohEff>0
                effect = effect+roiparams.cohoffset;
            end
        end
        effect = effect+inEff;

        cc_model(1,i,:) = data.canonical(data.cc.time(i)==data.utimes,:)*effect;

        res(i) = effect-data.cc.resp_(i);
    end
end

time_model = zeros(size(data.time.resp));

for i = 1:length(data.time.resp_)
    ccon = data.time.con(i);
    ccoh = data.time.coh(i);
    
    conEff = roiparams.sigmacon*(conModel(ccon,roiparams)-baseConResp);
    cohEff = roiparams.sigmacoh*(cohModel(ccoh,roiparams)-baseCohResp);
    
    if conEff==0 && cohEff==0 % no change! res=0
        res(length(data.time.resp_)+i) = 0;
    else
        if isfield(roiparams,'offset')
            effect = conEff+cohEff;
            if (conEff>0) || (cohEff>0)
                effect = effect+roiparams.offset;
            end
        else
            effect = conEff+cohEff;
            if conEff>0
                effect = effect+ roiparams.conoffset;
            end
            if cohEff>0
                effect = effect+roiparams.cohoffset;
            end
        end

        time_model(1,i,:) = data.canonical(data.time.time(i)==data.utimes,:)*effect;

        res(length(data.time.resp_)+i) = effect-data.time.resp_(i);
    end
end

if fixedParams.regularize
    res = [res 0.1*conModel(0:.1:1,roiparams) 0.1*cohModel(0:.1:1,roiparams)];
end
    
if 0%f>0
    figure(f)
    subplot(3,2,1:2);
    plot(res);
    subplot(3,2,3);
    x = 0:.01:1;
    plot(x,conModel(x,roiparams)*roiparams.sigmacon);
    subplot(3,2,4);
    plot(x,cohModel(x,roiparams)*roiparams.sigmacoh);
%     subplot(3,2,5);
%     hold on
%     plot(squeeze(fit.cc.resp)');
%     plot(squeeze(cc_model)');
    subplot(3,2,5:6);
    plot(fit.impulse);
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