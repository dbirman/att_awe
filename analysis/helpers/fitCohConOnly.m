function fitCohConOnly( nSIDs )

global fixedParams
%%
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
fixedParams.ROIs = rois;
rmode = {'resp','mresp'};
% subject | roi | model | index | value
cc = nan(length(nSIDs),length(rois),2,20,81);
for ni = 1:length(nSIDs)
    load(fullfile(datafolder,sprintf('s%04.0f_decon.mat',nSIDs(ni))));
    for ri = 1:length(rois)
        for rj = 1:length(rmode)
            cc(ni,ri,rj,:,:) = decondata.(rois{ri}).cc.(rmode{rj});
        end
    end
end

%% get indexes

conidx = decondata.V1.cc.conidxs;
cohidx = decondata.V1.cc.cohidxs;

cononly = cohidx==0;
cohonly = conidx==0.25;

con = cc(:,:,:,cononly,:);
coh = cc(:,:,:,cohonly,:);

fits = cell(size(nSIDs));
for ni = 1:length(nSIDs)
    fits{ni} = fitCCDecon(squeeze(con(ni,:,:,:,:)),squeeze(coh(ni,:,:,:,:)));
end

save(fullfile(datafolder,'avg_deconfits.mat'),'fits');

function fit = fitCCDecon(con,coh)

data.con = con;
data.coh = coh;
%% Format
% data is a roi * model * stimulus * resp matrix
% we will fit a single HRF, scaled by roi*stimulus

% Initialize parameters
% hrfparams.amp1 = 1;
% hrfparams.tau1 = [0.4 -inf inf];
% hrfparams.timelag1 = [1.5 0 3];
% hrfparams.amp2 = [-0.1 -inf 0];
% hrfparams.tau2 = [5 -inf inf];
% hrfparams.timelag2 = [-10 -inf inf];
% hrfparams.exponent = 8;
hrfparams = struct;
roiparams.con0 = 0;
roiparams.con1 = [1 0 inf];
roiparams.con2 = [2 0 inf];
roiparams.con3 = [3 0 inf];
roiparams.coh0 = 0;
roiparams.coh1 = [0.1 0 inf];
roiparams.coh2 = [0.2 0 inf];
roiparams.coh3 = [0.3 0 inf];
roiparams.coh4 = [0.4 0 inf];

global params

params.hrfparams = hrfparams;
params.roiparams = roiparams;

fit = fitModel(data);

function fit = fitModel(data)
%%
global fixedParams

[initparams,minparams,maxparams] = initParams;

f = figure;

load(fullfile(datafolder,'avg_hrf.mat'));
data.impulse = mean(hrfs);
% convolve the impulse and re-scale to 1
temp = conv([1 1 1 1 1],data.impulse);
temp = temp(1:length(data.impulse));
data.impulse = temp / max(temp);

optimParams = optimset('Algorithm','levenberg-marquardt','MaxIter',inf,'Display','off');
[bestparams, ~, ~, ~, ~, ~, ~] = lsqnonlin(@ccresid,initparams,minparams,maxparams,optimParams,data,fixedParams,f);

[~,fit] = ccresid(bestparams,data,fixedParams,-1);

fit.params = getParams(bestparams,fixedParams);
fit.roiparams = cell(size(fixedParams.ROIs));
fit.roiparams.con = zeros(length(fixedParams.ROIs),4);
fit.roiparams.coh = zeros(length(fixedParams.ROIs),5);
for ri = 1:length(fixedParams.ROIs)
    crp = getROIParams(fit.params,fixedParams.ROIs{ri});
    fit.roiparams.con(ri,:) = [crp.con0 crp.con1 crp.con2 crp.con3];
    fit.roiparams.coh(ri,:) = [crp.coh0 crp.coh1 crp.coh2 crp.coh3 crp.coh4];
end
fit.impulse = data.impulse;
% fit.impulse = cc_gamma(0.25:0.5:30.5,fit.params);

function [res,fit] = ccresid(params,data,fixedParams,f)
params = getParams(params,fixedParams);

t = 0.25:0.5:40.5;
% impulse = cc_gamma(t,params);
impulse = data.impulse;

model.con = zeros(size(data.con));
model.coh = zeros(size(data.coh));

res = zeros(1,length(t)*length(fixedParams.ROIs));
for ri = 1:length(fixedParams.ROIs)
    rp = getROIParams(params,fixedParams.ROIs{ri});
    
    con = squeeze(data.con(ri,1,:,:));
    coh = squeeze(data.coh(ri,1,:,:));
    
    dat = [con(1,:) con(2,:) con(3,:) con(4,:) coh(1,:) coh(2,:) coh(3,:) coh(4,:) coh(5,:)];
    lmodel = [impulse*rp.con0 impulse*rp.con1 impulse*rp.con2 impulse*rp.con3 impulse*rp.coh0 impulse*rp.coh1 impulse*rp.coh2 impulse*rp.coh3 impulse*rp.coh4];
    
    for i = 1:4
        model.con(ri,1,i,:) = impulse*rp.(sprintf('con%i',i-1));
    end
    for i = 1:5
        model.coh(ri,1,i,:) = impulse*rp.(sprintf('coh%i',i-1));
    end
    
    res((ri-1)*length(t)*9+1:(ri-1)*length(t)*9+length(t)*9) = dat-lmodel;
    
    if ri==1 && f>0
        figure(f) 
        clf
        hold on
        plot(lmodel,'b');
        plot(dat,'r');
    end
end

fit.data = data;
fit.model = model;
fit.t = t;
fit.tlength = length(t);

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