function fit = fitCCHRFModel( data , mode, pfit, dataopt, cvflag)
%CCROIMODEL Fit the contrast coherence model to an ROI
%
%   Dan Birman - Gardner Lab, Stanford University
%       May 10th, 2016
%
%   fit = fitCCHRFModel(data)
%
%   OUTPUT:

%% Setup

global fixedParams params
fixedParams = struct;
fixedParams.ROIs = data.ROIs;


%% Set data
if ~isempty(dataopt)
    data.cc.cresp = data.cc.(dataopt);
    data.time.cresp = data.time.(dataopt);
end

%% Permutation
if strfind(mode,'permute')
    % permute data
    perm = randperm(size(data.cc.cresp,2));
    data.cc.cresp = data.cc.cresp(:,perm,:);
    perm = randperm(size(data.time.cresp,2));
    data.time.cresp = data.time.cresp(:,perm,:);
end

%% Cross-validation
if cvflag
    disp('(fitCCHRFModel) CROSS-VALIDATION INITIATED');
    
    % split the data into train and test and run on each train, evaluated
    % on each test

    all = struct;
    all.y = []; all.y_ = [];
    all.train_r2 = zeros(39,8);
    all.test_r2 = zeros(39,8);
    
    ccn = length(data.cc.con);
    timen = length(data.time.con);
    total = ccn+timen;
    disppercent(-1/(total-1));
    % skip ci = 1, because that's the zero/zero condition
    for ci = 2:total
        % tdata will be used to train
        tdata = data;
        % test will be used to test
        test = data;
        % set the training data    
        if ci<=ccn
            % train
            train = setdiff(1:ccn,ci);
            % test is in the cc group
            tdata.cc.con = tdata.cc.con(train);
            tdata.cc.coh = tdata.cc.coh(train);
            tdata.cc.time = tdata.cc.time(train);
            tdata.cc.cresp = tdata.cc.cresp(:,train,:);
            % don't change the tdata time group
            % test
            test.cc.con = test.cc.con(ci);
            test.cc.coh = test.cc.coh(ci);
            test.cc.time = test.cc.time(ci);
            test.cc.cresp = test.cc.cresp(:,ci,:);
            % remove the test time group
            test.time.con = [];
            test.time.coh = [];
            test.time.time = [];
            test.time.cresp = test.time.cresp(:,[],:); % maintain size so the remaining code doesn't fail
        else
            % train
            train = setdiff(1:timen,ci-ccn);
            % test is in the time group
            tdata.time.con = tdata.time.con(train);
            tdata.time.coh = tdata.time.coh(train);
            tdata.time.time = tdata.time.time(train);
            tdata.time.cresp = tdata.time.cresp(:,train,:);
            % don't change the cc group
            % test
            test.time.con = test.time.con(ci-ccn);
            test.time.coh = test.time.coh(ci-ccn);
            test.time.time = test.time.time(ci-ccn);
            test.time.cresp = test.time.cresp(:,ci-ccn,:);
            % remove the test time group
            test.cc.con = [];
            test.cc.coh = [];
            test.cc.time = [];
            test.cc.cresp = test.cc.cresp(:,[],:); 
        end
        
        % you have to remove dataopt and set to [] otherwise this code
        % fails
        trainfit = fitCCHRFModel(tdata,mode,pfit,[],0); % no CV flag
        
        testfit = fitCCHRFModel(test,'predict',trainfit,[],0);
        
        all.y{ci-1} = testfit.y;
        all.y_{ci-1} = testfit.y_;
        
        all.train_r2(ci-1,:) = trainfit.r2;
        all.test_r2(ci-1,:) = testfit.r2;
        
        % clear explicitly
        clear tdata test train trainfit testfit
        disppercent((ci-1)/(total-1));
    end
    disppercent(inf);
    
    y = zeros(length(all.y)*81,8);
    y_ = y;
    for i = 1:length(all.y)
        for ri = 1:8
            y((i-1)*81+1:(i-1)*81+81,ri) = all.y{i}{ri};
            y_((i-1)*81+1:(i-1)*81+81,ri) = all.y_{i}{ri};
        end
    end
    all.y = y;
    all.y_ = y_;
    
    for ri = 1:8
        all.r2(ri) = myr2(y(:,ri),y_(:,ri));
    end

    % fit once to ALL available data
    fit = fitCCHRFModel(data,mode,pfit,dataopt,0);
    % merge the CV splits and compute r2
    fit.cv = all;
    return
end

%% If multiple ROIs, fit each individually
if length(data.ROIs)>1
    fit = struct;
    fit.y = {}; fit.y_ = {};
    for ri = 1:length(data.ROIs)
        disp(sprintf('(fitCCHRFModel) Running ROI: %s',data.ROIs{ri}));
        
        ndata = data;
        ndata.ROIs = ndata.ROIs(ri);
        ndata.cc.cresp = ndata.cc.cresp(ri,:,:);
        ndata.time.cresp = ndata.time.cresp(ri,:,:);
        if strfind(mode,'predict')
            ndata.params = pfit.roifit{ri}.params;
        end
        if strfind(mode,'fitexp')
            ndata.params = pfit.roifit{ri}.params;
        end
        % run model
        fit.roifit{ri} = fitCCHRFModel(ndata,mode,pfit,'',cvflag);
        % copy fit parameters
        fit.y{ri} = fit.roifit{ri}.y;
        fit.y_{ri} = fit.roifit{ri}.y_;
        fit.r2(ri) = fit.roifit{ri}.r2;
        fit.BIC(ri) = fit.roifit{ri}.BIC;
        fit.AIC(ri) = fit.roifit{ri}.AIC;
        fit.like(ri) = fit.roifit{ri}.like;
    end
   
    return
end

%% predict mode
if strfind(mode,'predict')
    hrfparams = struct;
    roiparams = struct;
    fixedParams.fitroi = 0;
    fixedParams.spkdec = 0;
    fixedParams.fitexp = 0;
    fixedParams.regularize = 0;
    if isfield(data.params,'cohalpha_0')
        fixedParams.intzero = 1;
    else
        fixedParams.intzero=0;
    end
    
    hrfparams.spkexp = 0;
    hrfparams.hrfexp = -0.451916488002049;
    
    fds = fields(data.params);
    for fi = 1:length(fds)
        roiparams.(fds{fi}) = data.params.(fds{fi});
    end
    

    params.hrfparams = hrfparams;
    params.roiparams = roiparams;
    
    fit = fitModel(data);
    
    return
end
%% parse mode:
hrfparams = struct;
roiparams = struct;
fixedParams.fitroi = 0;
fixedParams.spkdec = 0;
fixedParams.fitexp = 0;

if strfind(mode,'spkdec')
    % Spike rate decay
    hrfparams.spkexp = -0.451916488002049;%[-0.5 -inf 0];
    hrfparams.hrfexp = 0;
    fixedParams.spkdec = 1;
elseif strfind(mode,'spkhrfdec')
    hrfparams.spkexp = -0.451916488002049;
    hrfparams.hrfexp = -0.451916488002049;
else
    % HRF response decay
    hrfparams.spkexp = 0;
    hrfparams.hrfexp = -0.451916488002049; % adaptation exponent (for time)
end

fixedParams.numparams = 0;
if strfind(mode,'null')
    roiparams.conslope = 0;
    roiparams.conmodel = 1;
    roiparams.cohslope = 0;
    roiparams.cohmodel = 1;
    roiparams.nulloffset = [1 -inf inf];
elseif strfind(mode,'fitroi')
    if strfind(mode,'conlinear')
        roiparams.conslope = [1 -inf inf];
        roiparams.conmodel = 1;
    else
        roiparams.conRmax = [2 -inf inf];
        roiparams.conp = 0.3;%[ eps inf];
        roiparams.conq = 1.6;
        roiparams.conc50 = [0.5 eps 1-eps];
        roiparams.conmodel = 2;
    end
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
        roiparams.cohslope = [1 -inf inf];
        roiparams.cohmodel = 1;
    elseif strfind(mode,'naka')
        roiparams.cohRmax = [0.5 -inf inf];
        roiparams.cohp = 0.3;
        roiparams.cohq = 2;
        roiparams.cohc50 = [0.5 eps 1-eps];
        roiparams.cohmodel = 2;
    elseif strfind(mode,'exp')
        roiparams.cohalpha = [2 -inf inf];
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

if strfind(mode,'intzero')
    % Define two separate sets of paramters for the conditions with zero
    % change in one feature and change in the other feature and vice versa
    
    % Add the alpha/kappa parameters for the no contrast condition
    roiparams.cohalpha_0 = roiparams.cohalpha;
    roiparams.cohkappa_0 = roiparams.cohkappa;
    
    % Add the Rmax and c50 parameters for the no coherence condition
    roiparams.conRmax_0 = roiparams.conRmax;
    roiparams.conc50_0 = roiparams.conc50;
    
    fixedParams.intzero = 1;
else
    fixedParams.intzero = 0;
end

fixedParams.regularize=0;
if strfind(mode,'doreg')
    fixedParams.regularize = 1;
end

%% Parameter initialization

params.hrfparams = hrfparams;
params.roiparams = roiparams;

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
out_cc = zeros(1,size(data.cc.cresp,2));
for i=1:size(data.cc.cresp,2)
    curhrf = data.canonical(data.utimes==data.cc.time(i),:);
    out_cc(i) = curhrf'\squeeze(data.cc.cresp(1,i,:));
end
data.cc.cresp_ = out_cc;

out_time = zeros(1,size(data.time.cresp,2));
for i=1:size(data.time.time,2)
    curhrf = data.canonical(data.utimes==data.time.time(i),:);
    out_time(i) = curhrf'\squeeze(data.time.cresp(1,i,:));
end
data.time.cresp_ = out_time;

% f = figure;
f = 0;

if ~isempty(initparams)
    optimParams = optimset('Algorithm','trust-region-reflective','MaxIter',inf,'Display','off');
    [bestparams, ~, ~, ~, ~, ~, ~] = lsqnonlin(@hrfResidual,initparams,minparams,maxparams,optimParams,data,f,fixedParams);
else
    bestparams = initparams;
end

n = length(data.cc.cresp(:))+length(data.time.cresp(:));

[~,fit] = hrfResidual(bestparams,data,-1,fixedParams);
fit.SSE = sum(fit.rres.^2);
fit.BIC = n*log(fit.SSE/n) + length(bestparams)*log(n);
fit.AIC = n*log(fit.SSE/n) + length(bestparams)*2;
fit.like = log(fit.SSE/n);
fit.params = getParams(bestparams,fixedParams);
fit.params = getROIParams(fit.params,data.ROIs{1});
fit.ROIs = fixedParams.ROIs;

function [res, fit] = hrfResidual(params,data,f,fixedParams)

fit = struct;

% stimvol basecon lcon rcon basecoh lcoh rcoh timing task
if ~isstruct(params)
    params = getParams(params,fixedParams);
    roiparams = getROIParams(params,fixedParams.ROIs{1});
else
    roiparams = params;
end

fit.cc = data.cc;
fit.time = data.time;
fit.cc.model = zeros(size(fit.cc.cresp));
fit.time.model = zeros(size(fit.time.cresp));

res = zeros(1,length(data.cc.cresp_)+length(data.time.cresp_));

baseConResp = conModel(data.basecon,roiparams);
baseCohResp = cohModel(data.basecoh,roiparams);

if fixedParams.intzero
    roiparams_0 = roiparams;
    roiparams_0.conc50 = roiparams_0.conc50_0;
    roiparams_0.conRmax = roiparams_0.conRmax_0;
    roiparams_0.cohkappa = roiparams_0.cohkappa_0;
    roiparams_0.cohalpha = roiparams_0.cohalpha_0;
    baseConResp_0 = conModel(data.basecon,roiparams_0);
    baseCohResp_0 = cohModel(data.basecoh,roiparams_0);
end

cc_model = zeros(size(data.cc.cresp));

for i = 1:length(data.cc.cresp_)
    ccon = data.cc.con(i);
    ccoh = data.cc.coh(i);
    
    if fixedParams.intzero && (ccoh==0)
        conEff = conModel(ccon,roiparams_0)-baseConResp_0;
    else
        conEff = conModel(ccon,roiparams)-baseConResp;
    end
    if fixedParams.intzero && (ccon==0)
        cohEff = cohModel(ccoh,roiparams_0)-baseCohResp_0;
    else
        cohEff = cohModel(ccoh,roiparams)-baseCohResp;
    end
    inEff = roiparams.inbeta*conEff*cohEff;
    
    if conEff==0 && cohEff==0 && ~isfield(roiparams,'nulloffset') % no change! res=0
        res(i) = 0;
    else
        if isfield(roiparams,'nulloffset')
            effect = roiparams.nulloffset;
        elseif isfield(roiparams,'offset')
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

        res(i) = effect-data.cc.cresp_(i);
    end
end

time_model = zeros(size(data.time.cresp));

for i = 1:length(data.time.cresp_)
    ccon = data.time.con(i);
    ccoh = data.time.coh(i);
    
    if fixedParams.intzero && (ccoh==0)
        conEff = conModel(ccon,roiparams_0)-baseConResp_0;
    else
        conEff = conModel(ccon,roiparams)-baseConResp;
    end
    if fixedParams.intzero && (ccon==0)
        cohEff = cohModel(ccoh,roiparams_0)-baseCohResp_0;
    else
        cohEff = cohModel(ccoh,roiparams)-baseCohResp;
    end
    
    if conEff==0 && cohEff==0 && ~isfield(roiparams,'nulloffset') % no change! res=0
        res(length(data.time.cresp_)+i) = 0;
    else
        if isfield(roiparams,'nulloffset')
            effect = roiparams.nulloffset;
        elseif isfield(roiparams,'offset')
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

        res(length(data.time.cresp_)+i) = effect-data.time.cresp_(i);
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
    plot(x,conModel(x,roiparams));
    subplot(3,2,4);
    plot(x,cohModel(x,roiparams));
%     subplot(3,2,5);
%     hold on
%     plot(squeeze(fit.cc.cresp)');
%     plot(squeeze(cc_model)');
    subplot(3,2,5:6);
    plot(fit.impulse);
end

fit.cc.model = cc_model;
fit.time.model = time_model;
 
% Compute r^2

y = [data.cc.cresp(:); data.time.cresp(:)];
y_ = [cc_model(:); time_model(:)];

fit.y = y;
fit.y_ = y_;
fit.rres = y_-y;
fit.r2 = myr2(y,y_);
% fit.r = corrcoef([y y_]);
% fit.r2 = fit.r(1,2)^2;

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