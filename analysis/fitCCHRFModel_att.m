function fit = fitCCHRFModel_att( data , subj, mode, cvflag)
%CCROIMODEL Fit the contrast coherence model to an ROI
%
% Okay the new attention model works as follows: it loads the contrast /
% coherence response functions for this subject. Then it tries to test
% whether under the attention conditions the data is better fit by
% increasing the offset parameter (constant across trials) or a gain
% parameter (response gain effect), or both, or none of the above. 
%
if isempty(data)
    warning('No data available!!');
    fit = struct;
    return
end

%% Setup

global fixedParams params
fixedParams = struct;
fixedParams.ROIs = data.ROIs;
fixedParams.baseoffset = 0;

%% Load subject's existing model information
if ~isfield(data,'cc_fit')
    data.cc_fit = load(fullfile(datafolder,sprintf('s%04.0f_hrf.mat',subj)));
end

%% Cross-validation
if ~isstruct(cvflag) && cvflag>0
    disp('(fitCCHRFModel) CROSS-VALIDATION INITIATED');
    
    % split the data into train and test and run on each train, evaluated
    % on each test

    all = struct;
    all.y = []; all.y_ = [];
    all.train_r2 = zeros(1,8);
    all.test_r2 = zeros(1,8);
    
    total = 8;
    disppercent(-1/(total-1));
    
    idxs = find(~data.deltaidx);
    % skip ci = 1, because that's the zero/zero condition
    for ci = 1:total
        % tdata will be used to train
        tdata = data;
        % test will be used to test
        test = data;
        % set the training data  
        train = setdiff(1:total,ci);
        train = idxs(train);
        test_ = idxs(ci);
        
        tdata.conidx = tdata.conidx(train);
        tdata.cohidx = tdata.cohidx(train);
        tdata.taskidx = tdata.taskidx(train);
        tdata.deltaidx = tdata.deltaidx(train);
        tdata.conStim = tdata.conStim(train);
        tdata.cohStim = tdata.cohStim(train);
        tdata.betaCon = tdata.betaCon(:,:,train);
        tdata.betaCoh = tdata.betaCoh(:,:,train);
        
        test.conidx = test.conidx(test_);
        test.cohidx = test.cohidx(test_);
        test.taskidx = test.taskidx(test_);
        test.deltaidx = test.deltaidx(test_);
        test.conStim = test.conStim(test_);
        test.cohStim = test.cohStim(test_);
        test.betaCon = test.betaCon(:,:,test_);
        test.betaCoh = test.betaCoh(:,:,test_);
        
        % you have to remove dataopt and set to [] otherwise this code
        % fails
        trainfit = fitCCHRFModel_att(tdata,subj,mode,0); % no CV flag
        
        testfit = fitCCHRFModel_att(test,subj,mode,trainfit);
        
        all.y{ci} = testfit.y(:);
        all.y_{ci} = testfit.y_(:);
        
        % clear explicitly
        clear tdata test train trainfit testfit
        disppercent((ci-1)/(total-1));
    end
    disppercent(inf);
    
    y = zeros(8,16);
    y_ = y;
    for i = 1:length(all.y)
        y(i,:) = all.y{i};
        y_(i,:) = all.y_{i};
    end
    all.y = y(:);
    all.y_ = y_(:);
    
    for ri = 1:8
        all.r2(ri) = myr2(y(:,ri),y_(:,ri));
    end

    % fit once to ALL available data
    fit = fitCCHRFModel_att(data,subj,mode,0);
    % merge the CV splits and compute r2
    fit.cv = all;
    return
end

%% If multiple ROIs, fit each individually
if length(data.ROIs)>1
    for ri = 1:length(data.ROIs)
        disp(sprintf('(fitCCHRFModel) Running ROI: %s',data.ROIs{ri}));
        
        ndata = data;
        ndata.ROIs = ndata.ROIs(ri);
        ndata.cc_fit = ndata.cc_fit.cur.roifit{ri};
        ndata.betaCon = squeeze(ndata.betaCon(:,ri,:));
        ndata.betaCoh = squeeze(ndata.betaCoh(:,ri,:));
        if isstruct(cvflag)
            fit.roifit{ri} = fitCCHRFModel_att(ndata,subj,mode,cvflag.roifit{ri});
        else
            fit.roifit{ri} = fitCCHRFModel_att(ndata,subj,mode,0);
        end
        fit.y(ri,:) = fit.roifit{ri}.y;
        fit.y_(ri,:) = fit.roifit{ri}.y_;
        fit.r2(ri) = fit.roifit{ri}.r2;
        fit.AIC(ri) = fit.roifit{ri}.AIC;
    end
   
    return
end

%%
if isstruct(cvflag)
    % evaluate test model and return prediction
    fixedParams.x = 0:.001:1;
    fixedParams.con = conModel(fixedParams.x,data.cc_fit.params);
    fixedParams.coh = cohModel(fixedParams.x,data.cc_fit.params);
    p = cvflag.roiparams;
    hrfparams = struct; roiparams = struct;
    if isfield(p,'offset_con')
        fixedParams.offset = 2;
        roiparams.offset_con = p.offset_con;
        roiparams.offset_coh = p.offset_coh;
    elseif p.offset_shift > 0
        fixedParams.offset = 1;
        roiparams.offset_shift = p.offset_shift;
    else
        fixedParams.offset = 0;
        roiparams.offset_shift = p.offset_shift;
    end

    roiparams.con_congain = p.con_congain;
    roiparams.con_cohgain = p.con_cohgain;
    roiparams.coh_congain = p.coh_congain;
    roiparams.coh_cohgain = p.coh_cohgain;
    
    roiparams.conmodel = 4;
    roiparams.cohmodel = 4;
    params.hrfparams = hrfparams;
    params.roiparams = roiparams;
    
    fit = fitModel(data);
    return
end

%%

fixedParams.x = 0:.001:1;
fixedParams.con = conModel(fixedParams.x,data.cc_fit.params);
fixedParams.coh = cohModel(fixedParams.x,data.cc_fit.params);

%% parse mode:
% "fithrf" - just fit the HRF using all trials (but no effects)
% "fitroi" - use a computed HRF to fit the trials
% functions
% "fitall" - run fithrf, fitroi, and fitatt and return the full fit
hrfparams = struct;
roiparams = struct;

hrfparams.baseoffset = data.cc_fit.params.offset;

fixedParams.offset = 0;

if ~isempty(strfind(mode,'doublebaseline'))
    roiparams.offset_con = [0 -inf inf];
    roiparams.offset_coh = [0 -inf inf];
    fixedParams.offset = 2;
elseif ~isempty(strfind(mode,'nobaseline'))
    roiparams.offset_shift = 0;
else
    roiparams.offset_shift = [0 -inf inf];
    fixedParams.offset = 1;
end

if ~isempty(strfind(mode,'nogain'))
    roiparams.con_congain = 1;
    roiparams.con_cohgain = 1;
    roiparams.coh_congain = 1;
    roiparams.coh_cohgain = 1;
    roiparams.conmodel = 4;
    roiparams.cohmodel = 4;
else
    roiparams.con_congain = [1 -inf inf];
    roiparams.con_cohgain = [1 -inf inf];
    roiparams.coh_congain = [1 -inf inf];
    roiparams.coh_cohgain = [1 -inf inf];
    roiparams.conmodel = 4;
    roiparams.cohmodel = 4;
end

%% Parameter initialization

params.hrfparams = hrfparams;
params.roiparams = roiparams;

%% fit
fit = fitModel(data);

fit.mode = mode;

function fit = fitModel(data)
%%
% Fit to the mean timeseries from the top R^2 values
global fixedParams
[initparams,minparams,maxparams] = initParams;

f = -inf;

if length(initparams)>0
    optimParams = optimset('Algorithm','levenberg-marquardt','MaxIter',inf,'Display','off');
    [bestparams, ~, ~, ~, ~, ~, ~] = lsqnonlin(@hrfResidual,initparams,minparams,maxparams,optimParams,data);
else
    bestparams=initparams;
end

[res,fit] = hrfResidual(bestparams,data);
fit.SSE = sum(res.^2);
RSS = fit.SSE;
n = sum(data.deltaidx==0);
fit.AIC = n*log(RSS/n) + length(bestparams)*2;
fit.params = getParams(bestparams,fixedParams);
fit.roiparams = getROIParams(fit.params,data.ROIs{1});
fit.ROIs = fixedParams.ROIs;

% draw response functions
fit.conresp = fixedParams.con + fixedParams.baseoffset;
if fixedParams.offset==2
    fit.conresp_con = fixedParams.con*fit.roiparams.con_congain + fixedParams.baseoffset + fit.roiparams.offset_con;
    fit.conresp_coh = fixedParams.con*fit.roiparams.con_cohgain + fixedParams.baseoffset + fit.roiparams.offset_con;
else
    fit.conresp_con = fixedParams.con*fit.roiparams.con_congain + fixedParams.baseoffset + fit.roiparams.offset_shift;
    fit.conresp_coh = fixedParams.con*fit.roiparams.con_cohgain + fixedParams.baseoffset + fit.roiparams.offset_shift;
end
fit.condata_con = data.betaCon(logical((data.taskidx==2).*(data.deltaidx==0)));
fit.condata_coh = data.betaCon(logical((data.taskidx==1).*(data.deltaidx==0)));

fit.cohresp = fixedParams.coh + fixedParams.baseoffset;
if fixedParams.offset==2
    fit.cohresp_con = fixedParams.coh*fit.roiparams.coh_congain + fixedParams.baseoffset + fit.roiparams.offset_coh;
    fit.cohresp_coh = fixedParams.coh*fit.roiparams.coh_cohgain + fixedParams.baseoffset + fit.roiparams.offset_coh;
else
    fit.cohresp_con = fixedParams.coh*fit.roiparams.coh_congain + fixedParams.baseoffset + fit.roiparams.offset_shift;
    fit.cohresp_coh = fixedParams.coh*fit.roiparams.coh_cohgain + fixedParams.baseoffset + fit.roiparams.offset_shift;
end
fit.cohdata_con = data.betaCoh(logical((data.taskidx==2).*(data.deltaidx==0)));
fit.cohdata_coh = data.betaCoh(logical((data.taskidx==1).*(data.deltaidx==0)));


%%%%%%%%%%%%
%%% HRF RESIDUAL
%%%%%%%%%%%%%%
function [res, fit] = hrfResidual(params,data)
%%
global fixedParams
fit = struct;
params = getParams(params,fixedParams);

%% Compute
roiparams = getROIParams(params,data.ROIs{1});

% compute gains
con_conresp = fixedParams.con * roiparams.con_congain;
con_cohresp = fixedParams.con * roiparams.con_cohgain;
% coh
coh_conresp = fixedParams.coh * roiparams.coh_congain;
coh_cohresp = fixedParams.coh * roiparams.coh_cohgain;

% compute residual
res = zeros(1,2*length(data.conidx));
y_ = zeros(1,2*length(data.conidx));
for i = 1:length(data.conidx)
    
    % contrast!
    pos = find(fixedParams.x>=data.conidx(i),1);
    if data.taskidx(i)==1
        resp = con_cohresp(pos);
    else
        resp = con_conresp(pos);
    end
    if fixedParams.offset==2
        if data.taskidx(i)==1
            resp = resp + fixedParams.baseoffset + roiparams.offset_coh;
        else
            resp = resp + fixedParams.baseoffset + roiparams.offset_con;
        end
    else
        resp = resp + fixedParams.baseoffset + roiparams.offset_shift;
    end
    
    res(i) = resp - data.betaCon(i);
    y_(i) = resp;
    clear resp
    
    % coherence!
    pos = find(fixedParams.x>=data.cohidx(i),1);
    if data.taskidx(i)==1
        resp = coh_cohresp(pos);
    else
        resp = coh_conresp(pos);
    end
    if fixedParams.offset==2
        if data.taskidx(i)==1
            resp = resp + fixedParams.baseoffset + roiparams.offset_coh;
        else
            resp = resp + fixedParams.baseoffset + roiparams.offset_con;
        end
    else
        resp = resp + fixedParams.baseoffset + roiparams.offset_shift;
    end
    
    res(length(data.conidx)+i) = resp - data.betaCoh(i);
    y_(length(data.conidx)+i) = resp;
end

% we only compute for the non-delta (8 values)
res = res(logical([data.deltaidx==0 data.deltaidx==0]));

fit.y = [data.betaCon; data.betaCoh];
fit.y_ = y_;

fit.r2 = myr2(fit.y,fit.y_);

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