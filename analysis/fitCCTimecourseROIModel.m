function fit = fitCCTimecourseROIModel( data , mode, fit)
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
%   Where Canonical is a cc_gamma function with three parameters: exponent,
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

%% Choose tSeries
if strfind(mode,'useprf')
    disp('(roimodel) Using pRF based averages');
    tSeriesname = 'tSeries';
else
    disp('(roimodel) Using top 25 voxel averages (per ROI)');
    tSeriesname = 'rtSeries25';
end

%% Concatenate Sessions if Necessary

if iscell(data)
    disp('(roimodel) Concatenating what appear to be different sessions...');
    % to run the model we really only need:
    % data.tSeries
    % data.design
    % data.runtrans
    % but we need these to be concatenated across the different sessions.
    % We have to be a little careful here that we make sure we add the
    % lengths of each run correctly or we'll screw ourselves over.
    data_old = data;
    data = struct;
    data.ROIs = data_old{1}.ROIs;
    % (we can do this the slow way, perf doesn't really matter
    data.tSeries = cell(1,length(data_old{1}.ROIs));
    length_sofar = 0;
    for ri = 1:length(data.tSeries), data.tSeries{ri} = []; end
    data.design = [];
    data.runtrans = [];
    for si = 1:length(data_old)
        for ri = 1:length(data_old{si}.tSeries)
            data.tSeries{ri} = [data.tSeries{ri} data_old{si}.(tSeriesname){ri}];
        end
        % tweak the SV by adding 
        cdes = data_old{si}.design;
        cdes(:,1) = cdes(:,1) + length_sofar;
        data.design = [data.design ; cdes];
        data.runtrans = [data.runtrans ; data_old{si}.runtrans+length_sofar];
        length_sofar = length_sofar + length(data_old{si}.tSeries{1});
    end
end

%% Drop timing
if ~isempty(strfind(mode,'droptiming'))
    disp('%% All Timing Responses Dropped %%');
    idxs = data.design(:,8)==5;
    data.design = data.design(idxs,:);
end

%% Cross-Validation
% if strfind(mode,'crossval')
%     disp('(roimodel) Cross-Validation (per run) is now running.');
%         
%     for i = 1:size(data.runtrans,1)
%         disp('(roimodel) Setting aside run %i as the test set',i);
%         testdata = struct; traindata=struct;
%         testdata.runtrans = data.runtrans(i,:);
%         others = 1:size(data.runtrans,1);
%         others = others(others~=i);
%         traindata.runtrans = data.runtrans(others,:);
%     end
%     stop = 1;
% end

%% Setup

if ~iscell(data.tSeries)
    data.tSeries{1} = data.tSeries;
end

global fixedParams
fixedParams = struct;
fixedParams.disp = 1;
fixedParams.diff = 1;
fixedParams.ROIs = data.ROIs;
fixedParams.fitting = 0;

if ~isempty(strfind(mode,'fitatt'))
    disp('Fitting by attention condition:');
    if ~isempty(strfind(mode,'fitatt=1'))
        disp('Motion');
        data.design = data.design(data.design(:,9)==1,:);
    else
        disp('Contrast');
        data.design = data.design(data.design(:,9)==2,:);
    end
end

%% If fitroi or fitatt, run by ROI
if ~isempty(strfind(mode,'fitroi'))
    fixedParams.fitting = 1;
    nfit = fit;
    if length(data.ROIs)>1
        nfit.r2 = zeros(size(data.ROIs));
        nfit.likelihoods = zeros(size(data.ROIs));
        nfit.chi2s = zeros(size(data.ROIs));
        for ri = 1:length(data.ROIs)
            disp(sprintf('(roimodel) Fitting %s',data.ROIs{ri}));
            rdata = data;
            rdata.ROIs = rdata.ROIs(ri);
            rdata.tSeries = rdata.tSeries(ri);
            rfit = fit;
%             rfit.ROIs = rfit.ROIs(ri);
            rfit.roiparams = rfit.roiparams(ri);
%             rfit.model = rfit.model(ri);
%             rfit.tSeries = rfit.tSeries(ri);
            rfit = fitCCTimecourseROIModel(rdata,mode,rfit);
            nfit.r2(ri) = rfit.r2;
            nfit.model{ri} = rfit.model{1};
            nfit.roiparams{ri} = rfit.roiparams{1};
            nfit.likelihoods(ri) = rfit.likelihood;
%             nfit.chi2s(ri) = rfit.chi2;
            afields = fields(rfit);
            for ai = 1:length(afields)
                if strfind(afields{ai},data.ROIs{ri})
                    % copy field
                    nfit.(afields{ai}) = rfit.(afields{ai});
                end
            end
        end
        fit = nfit;
        return
    end
end

%% parse mode:
% "fithrf" - just fit the HRF using all trials (but no effects)
% "fitroi" - use a computed HRF to fit the trials
% functions
% "fitall" - run fithrf, fitroi, and fitatt and return the full fit
hrfparams = struct;
roiparams = struct;
fixedParams.fithrf = 0;
fixedParams.fitroi = 0;
fixedParams.fitatt = 0;
fixedParams.refit = 0;
if strfind(mode,'refit')
    disp('Refitting HRF parameters');
    % add the min/max to the hrf params so that they will be fit
    hrfparams = copyhrfparams_refit(fit);
    % copy all the roi parameters
    roiparams = copyroiparams(fit);
    % copy all of these into hrf params, and make roiparams empty
    f = fields(roiparams);
    for fi = 1:length(f)
        hrfparams.(f{fi}) = roiparams.(f{fi});
    end
    roiparams = struct;
    fixedParams.refit = 1;
elseif strfind(mode,'fithrf')
    % HRF Parameters
    hrfparams.amp1 = 1;
    hrfparams.tau1 = [0.4 -inf inf];
    hrfparams.timelag1 = [1.5 0 3];
    hrfparams.amp2 = [-eps -inf 0];
    hrfparams.tau2 = [0.4 -inf inf];
    hrfparams.timelag2 = [4 0 9];
    hrfparams.exponent = 7;
    hrfparams.offset = [0 0 inf];
    % estimates how much the hrf drops off per stimvol
    hrfparams.adaptation = [1]; % NOT USING CURRENTLY
    roiparams.betas = [0.5 0 inf];
    fixedParams.fithrf = 1;
elseif strfind(mode,'fitroi')
    % get hrf params
    hrfparams = copyhrfparams(fit);
    % Offset
    if isfield(hrfparams,'offset')
        hrfparams = rmfield(hrfparams,'offset');
    end
    roiparams.offset = [0 -inf inf];
    % Contrast Function Parameters
    roiparams.conalpha = [1 -inf inf];
    roiparams.conkappa = [1 -inf inf];
    roiparams.conmodel = 3;
    % Coherence Function Parameters
    roiparams.cohalpha = [1 -inf inf];
    roiparams.cohkappa = [1 -inf inf];
    roiparams.cohmodel = 3;
    % run type
    fixedParams.fitroi = 1;
end

%% Parameter initialization
global params

params.hrfparams = hrfparams;
params.roiparams = roiparams;

fixedParams.sstot = 0;
for ti = 1:length(data.tSeries)
    data.tSeries{ti} = (data.tSeries{ti}-1)*100; % move into zero mean and 1% space, easier for interpretation later
    fixedParams.sstot = fixedParams.sstot + sum((data.tSeries{ti}-repmat(mean(data.tSeries{ti},2),1,size(data.tSeries{ti},2))).^2,2);
end

%% setup plots
if fixedParams.fithrf
    concatInfo.runTransition = data.runtrans;
    fixedParams.concatInfo = concatInfo;
    fixedParams.sv = {data.design(:,1)};
end

%% build timeseries mask
sv_view = 0.5:0.5:20;
tmin = 3;
tmax = 9;
sv_view = find(logical(sv_view>=tmin).*logical(sv_view<=tmax));
mask = zeros(size(data.tSeries{1}));
for run = 1:size(data.runtrans,1)
    cdesign = fil(data.design,1,'>=',data.runtrans(run,1));
    cdesign = fil(cdesign,1,'<=',data.runtrans(run,2));

    for si = 1:size(cdesign,1)
        % okay, for each stimvol, place its effect
        sv = cdesign(si,1);
        csv_view = min(sv+sv_view,data.runtrans(run,2));
        if ~all(csv_view==data.runtrans(run,2))
            mask(csv_view) = 1;
        end
    end
end
if strfind(mode,'nomask')
    mask = ones(size(data.tSeries{1}));
end
disp('Ignoring mask');
mask = ones(size(data.tSeries{1}));
fixedParams.mask = logical(mask);    

%% fit HRF
fit = fitModel(data);

fit.tSeriesname = tSeriesname;
fit.mode = mode;
fit.design = data.design;
fit.runtrans = data.runtrans;
if isfield(data,'basecon')
    fit.basecon = data.basecon;
    fit.basecoh = data.basecoh;
end

% reset fixedParams so it doesn't fuck us over later
fixedParams.fitting = 0;

function hrfparams = copyhrfparams(fit)
hrfparams.amp1 = fit.params.amp1;
if ~(hrfparams.amp1==1), warning('Amplitude is wrong?');keyboard; end
hrfparams.tau1 = fit.params.tau1;
hrfparams.timelag1 = fit.params.timelag1;
hrfparams.amp2 = fit.params.amp2;
hrfparams.tau2 = fit.params.tau2;
hrfparams.timelag2 = fit.params.timelag2;
hrfparams.exponent = fit.params.exponent;
hrfparams.adaptation = fit.params.adaptation;

function hrfparams = copyhrfparams_refit(fit)
hrfparams.amp1 = fit.params.amp1; % should be one
if ~(hrfparams.amp1==1), warning('Amplitude is wrong?');keyboard; end
hrfparams.tau1 = [fit.params.tau1 -inf inf];
hrfparams.timelag1 = [fit.params.timelag1 0 3];
hrfparams.amp2 = [fit.params.amp2 -inf 0];
hrfparams.tau2 = [fit.params.tau2 -inf inf];
hrfparams.timelag2 = [fit.params.timelag2 0 6];
hrfparams.exponent = fit.params.exponent;
hrfparams.adaptation = fit.params.adaptation;
    
function roiparams = copyroiparams(fit)
global fixedParams
% generic info
p = fit.roiparams{1};
f = fields(p);
roiparams = struct;
% copy all, or copy only one set from p
if length(fit.roiparams)>1
    % roiparams is a cell, so copy all rois separately (and rename them)
    for ri = 1:length(fit.ROIs)
        localparams = fit.roiparams{ri};
        for fi = 1:length(f)
            cfield = f{fi};
%             roiparams.(sprintf('%s%s',fit.ROIs{ri},cfield)) = [localparams.(cfield) -inf inf];
            roiparams.(sprintf('%s%s',fit.ROIs{ri},cfield)) = [localparams.(cfield)];
        end
        roiparams.(sprintf('%sconmodel',fit.ROIs{ri})) = localparams.conmodel;
        roiparams.(sprintf('%scohmodel',fit.ROIs{ri})) = localparams.cohmodel;
    end
else
    % only one set
    for fi = 1:length(f)
        cfield = f{fi};
        if fixedParams.refitroi
            roiparams.(cfield) = [p.(cfield) -inf inf];
        else
            roiparams.(cfield) = p.(cfield);
        end
    end
    % overwrite so that these don't get fit
    roiparams.conmodel = p.conmodel;
    roiparams.cohmodel = p.cohmodel;
end

function fit = fitModel(data)

% Fit to the mean timeseries from the top R^2 values
global fixedParams
[initparams,minparams,maxparams] = initParams;

if fixedParams.disp
    f = figure;
else
    f = -inf;
end

optimParams = optimset('Algorithm','trust-region-reflective','MaxIter',inf,'Display','off');
[bestparams, ~, ~, ~, ~, ~, ~] = lsqnonlin(@hrfResidual,initparams,minparams,maxparams,optimParams,data.tSeries,data.design,data.runtrans,f,fixedParams);

[res,fit] = hrfResidual(bestparams,data.tSeries,data.design,data.runtrans,0,fixedParams);
fit.SSE = sum(res.^2);
% fit.BIC = n*log(RSS/n) + params(gi)*log(n);
fit.params = getParams(bestparams,fixedParams);
fit.roiparams = cell(size(fixedParams.ROIs));
for ri = 1:length(fixedParams.ROIs)
    fit.roiparams{ri} = getROIParams(fit.params,fixedParams.ROIs{ri});
end
fit.ROIs = fixedParams.ROIs;

function [res, fit] = hrfResidual(params,tSeries,design,runtrans,f,fixedParams)

fit = struct;

% stimvol basecon lcon rcon basecoh lcoh rcoh timing task
params = getParams(params,fixedParams);

t = 0.25:0.5:50.5;
impulse = cc_gamma(t,params);
fit.impulse = impulse;
fit.t = t;

fit.model = cell(size(fixedParams.ROIs));
fit.tSeries = tSeries;
timepoints = length(tSeries{1});
res = zeros(1,timepoints*length(fixedParams.ROIs));
localres = zeros(1,timepoints*length(fixedParams.ROIs));
for ri = 1:length(fixedParams.ROIs)
    roiparams = getROIParams(params,fixedParams.ROIs{ri});
    % pick the indexes to use
    if strcmp(fixedParams.ROIs{ri}(1),'l')
        conidx = 4; % use right
        cohidx = 7;
    else
        conidx = 3;
        cohidx = 6;
    end
    ctSeries = tSeries{ri};
    roimodel = zeros(size(ctSeries));
    % get just this ROIs parameters and strip out the ROI name itself
    for run = 1:size(runtrans,1)
        cdesign = fil(design,1,'>=',runtrans(run,1));
        cdesign = fil(cdesign,1,'<=',runtrans(run,2));

        for si = 1:size(cdesign,1)
            % okay, for each stimvol, place its effect
            sv = cdesign(si,1);
            if ((cdesign(si,conidx)-cdesign(si,2))+(cdesign(si,cohidx)-cdesign(si,5)))>0
                if fixedParams.fitroi || fixedParams.fitatt || fixedParams.refit
                    if fixedParams.fitatt
                        coneff = conModel(cdesign(si,conidx)-cdesign(si,2),roiparams,cdesign(si,9),1);
                        coheff = cohModel(cdesign(si,cohidx)-cdesign(si,5),roiparams,cdesign(si,9),1);
                    else
                        coneff = conModel(cdesign(si,conidx)-cdesign(si,2),roiparams,0,1);
                        coheff = cohModel(cdesign(si,cohidx)-cdesign(si,5),roiparams,0,1);
                    end
                    effect = coneff+coheff;
                 
                    % adjust timing
                    if cdesign(si,8)>=1
                        idxs = sv:min(runtrans(run,2),sv+(cdesign(si,8)-1));
                    else
                        idxs = sv;
                        effect = effect * cdesign(si,8);
                    end
                    
%                     afunc = ones(size(idxs));
%                     roimodel(:,idxs) = roimodel(:,idxs) + effect + roiparams.offset;
                    roimodel(:,idxs) = roimodel(:,idxs) + effect;
                    roimodel(:,idxs(1)) = roimodel(:,idxs(1)) + roiparams.offset;

                else
                    % we're just fitting the HRF model
                    if cdesign(si,8)>=1
                        idxs = sv:min(runtrans(run,2),sv+(cdesign(si,8)-1));
                        effect = roiparams.betas;
                    else
                        idxs = sv;
                        effect = roiparams.betas*cdesign(si,8);
                    end
                    roimodel(:,idxs) = roimodel(:,idxs) + effect; % don't use offset, correlated to betas
                    roimodel(:,idxs(1)) = roimodel(:,idxs(1)) + params.offset;
                end
            end
        end
    end
    cmt = conv(impulse,roimodel);
    fit.model{ri} = cmt(1:size(ctSeries,2));
    res((ri-1)*timepoints+1:(ri-1)*timepoints+timepoints) = fixedParams.mask .* (ctSeries-fit.model{ri});
    % save a local copy of the residual just to use for the R^2 calculation
    localres((ri-1)*timepoints+1:(ri-1)*timepoints+timepoints) = (ctSeries-fit.model{ri});
end

% this isn't used by lsqnonlin so it's safe to use the non-masked versions
ssres = sum(localres.^2);
fit.r2 = 1 - ssres/fixedParams.sstot;

fit.likelihood = ssres;

if false % f>0
    figure(f)
    plot(impulse);
    return
    
    if fixedParams.refithrf
        figure(f)
        plot(impulse);
%         disp(sprintf('Skipping plot: ss = %4.2f',fit.likelihood));
%         disp(sprintf('Adaptation: %0.2f',params.adaptation));
    elseif fixedParams.fithrf
%         disp(sprintf('Skipping plot: ss = %4.2f',fit.likelihood));
%         disp(sprintf('Adaptation: %0.2f Offset: %1.2f',params.adaptation,params.offset));
        return
        figure(f)
        clf(f), hold on
        curd = constructD(fit.tSeries{1}/100+1,fixedParams.sv,0.5,50,fixedParams.concatInfo,'none','deconv',0);
        decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
        decon = rmfield(decon,'scm');
        decon = rmfield(decon,'covar');
        plot(decon.ehdr,'o');
        curd = constructD(fit.model{1}/100+1,fixedParams.sv,0.5,50,fixedParams.concatInfo,'none','deconv',0);
        decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
        decon = rmfield(decon,'scm');
        decon = rmfield(decon,'covar');
        plot(decon.ehdr);
        title(sprintf('Amp1: %1.2f Tau1: %1.2f TL1: %1.2f Amp2: %1.2f Tau2: %1.2f TL2: %1.2f ',params.amp1,params.tau1,params.timelag1,params.amp2,params.tau2,params.timelag2));
    else
        figure(f)
        subplot(211)
        hold on
        plot(ctSeries(1:1000),'b');
        plot(fit.model{end}(1:1000),'r');
        title(sprintf('R^2: %0.2f',fit.r2));
        subplot(212), hold on
        cmap = brewermap(7,'PuOr');
        x = 0:.01:1;
        plot(x,conModel(x,roiparams),'Color',cmap(2,:));
        plot(x,cohModel(x,roiparams),'Color',cmap(6,:));
    end
end

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