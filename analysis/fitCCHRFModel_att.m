function fit = fitCCHRFModel_att( data , subj, mode)
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

global fixedParams
fixedParams = struct;
fixedParams.ROIs = data.ROIs;

%% Load subject's existing model information
if ~isfield(data,'cc_fit')
    data.cc_fit = load(fullfile(datafolder,sprintf('s%04.0f_hrf.mat',subj)));
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
        
        fit.roifit{ri} = fitCCHRFModel_att(ndata,subj,mode);
        fit.r2(ri) = fit.roifit{ri}.r2;
        fit.BIC(ri) = fit.roifit{ri}.BIC;
    end
   
    return
end

%%

fixedParams.x = 0:.01:1;
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
global params

params.hrfparams = hrfparams;
params.roiparams = roiparams;

fixedParams.sstot = sum([data.betaCon(data.deltaidx==0)' data.betaCoh(data.deltaidx==0)'].^2);  

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
fit.BIC = n*log(RSS/n) + length(bestparams)*log(n);
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
for i = 1:length(data.conidx)
    
    % contrast!
    pos = find(fixedParams.x>=data.conidx(i),1);
    if data.taskidx(i)==1
        resp = con_cohresp(pos);
    else
        resp = con_conresp(pos);
    end
    if fixedParams.offset==2
        resp = resp + fixedParams.baseoffset + roiparams.offset_con;
    else
        resp = resp + fixedParams.baseoffset + roiparams.offset_shift;
    end
    
    res(i) = resp - data.betaCon(i);
    clear resp
    
    % coherence!
    pos = find(fixedParams.x>=data.cohidx(i),1);
    if data.taskidx(i)==1
        resp = coh_cohresp(pos);
    else
        resp = coh_conresp(pos);
    end
    if fixedParams.offset==2
        resp = resp + fixedParams.baseoffset + roiparams.offset_coh;
    else
        resp = resp + fixedParams.baseoffset + roiparams.offset_shift;
    end
    
    res(length(data.conidx)+i) = resp - data.betaCoh(i);
end

% we only compute for the non-delta (8 values)
res = res(logical([data.deltaidx==0 data.deltaidx==0]));

%% Finalize

% this isn't used by lsqnonlin so it's safe to use the non-masked versions
ssres = sum(res.^2);
fit.r2 = 1 - ssres/fixedParams.sstot;

fit.likelihood = ssres;

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