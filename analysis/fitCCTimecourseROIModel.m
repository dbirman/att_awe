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
    if strfind(mode,'useprf')
        disp('(roimodel) Using pRF based averages');
    else
        disp('(roimodel) Using top 25 voxel averages (per ROI)');
    end
    for si = 1:length(data_old)
        for ri = 1:length(data_old{si}.tSeries)
            if strfind(mode,'useprf')
                data.tSeries{ri} = [data.tSeries{ri} data_old{si}.tSeries{ri}];
            else
                data.tSeries{ri} = [data.tSeries{ri} data_old{si}.rtSeries{ri}];
            end
        end
        % tweak the SV by adding 
        cdes = data_old{si}.design;
        cdes(:,1) = cdes(:,1) + length_sofar;
        data.design = [data.design ; cdes];
        data.runtrans = [data.runtrans ; data_old{si}.runtrans+length_sofar];
        length_sofar = length_sofar + length(data_old{si}.tSeries{1});
    end
end

%% Setup

if ~iscell(data.tSeries)
    data.tSeries{1} = data.tSeries;
end

global fixedParams
fixedParams = struct;
fixedParams.disp = 1;
fixedParams.diff = 1;
fixedParams.ROIs = data.ROIs;

%% If fitroi or fitatt, run by ROI
if ~isempty(strfind(mode,'fitroi'))||~isempty(strfind(mode,'fitatt'))
    if length(data.ROIs)>1
        fit.r2 = zeros(size(data.ROIs));
        for ri = 1:length(data.ROIs)
            disp(sprintf('(roimodel) Fitting %s',data.ROIs{ri}));
            rdata = data;
            rdata.ROIs = rdata.ROIs(ri);
            rdata.tSeries = rdata.tSeries(ri);
            rfit = fitCCTimecourseROIModel(rdata,mode,fit);
            fit.r2(ri) = rfit.r2;
            fit.model{ri} = rfit.model{1};
            fit.roiparams{ri} = rfit.roiparams{1};
            afields = fields(rfit);
            for ai = 1:length(afields)
                if strfind(afields{ai},data.ROIs{ri})
                    % copy field
                    fit.(afields{ai}) = rfit.(afields{ai});
                end
            end
        end
        return
    end
end

%% parse mode:
% "fithrf" - just fit the HRF using all trials (but no effects)
% "fitroi" - use a computed HRF to fit the trials
% "fitatt" - use the computed HRF and ROI data to fit the attention
% functions
% "fitall" - run fithrf, fitroi, and fitatt and return the full fit
hrfparams = struct;
roiparams = struct;
fixedParams.fithrf = 0;
fixedParams.fitroi = 0;
fixedParams.fitatt = 0;
if strfind(mode,'fitall')
    warning('not implemented');
    keyboard
elseif strfind(mode,'fithrf')
    % HRF Parameters
    hrfparams.amp1 = 1;
    hrfparams.tau1 = [0.45 -inf inf];
    hrfparams.timelag1 = [.7 0 3];
    hrfparams.amp2 = [-0.25 -inf 0];
    hrfparams.tau2 = [1.7 -inf inf];
    hrfparams.timelag2 = [0 0 6];
    hrfparams.exponent = 7;
    hrfparams.n = 1;
    hrfparams.offset = 0;
    roiparams.betas = [1 0 inf];
    fixedParams.fithrf = 1;
elseif strfind(mode,'fitroi')
    % get hrf params
    hrfparams = copyhrfparams(fit);
    % add a global offset
    hrfparams.offset = [.05 -inf inf];
    % Contrast Function Parameters
    roiparams.n = 1;
    roiparams.Rmax = [.1 0 inf];
    roiparams.c50 = [0.5 0 1];
    % Coherence Function Parameters
    roiparams.slope = [0.1 -inf inf]; 
    fixedParams.fitroi = 1;
elseif strfind(mode,'fitatt')  
    hrfparams = copyhrfparams(fit);
    roiparams = copyroiparams(fit);
    if strfind(mode,'gain')
        roiparams.conatt = [1 0 inf];
        roiparams.conunatt = [1 0 inf];
        roiparams.cohatt = [1 0 inf];
        roiparams.cohunatt = [1 0 inf];
        fixedParams.att = 1;
    elseif strfind(mode,'add')
        roiparams.conatt = [0 -inf inf];
        roiparams.cohatt = [0 -inf inf];
        roiparams.conunatt = [0 -inf inf];
        roiparams.cohunatt = [0 -inf inf];
        fixedParams.att = 0;
    end
    if isempty(mode) || ~isempty(strfind(mode,'noatt'))
        disp('(roimodel) Not fitting attention');
        roiparams.conatt = 1;
        roiparams.cohatt = 1;
        roiparams.conunatt = 1;
        roiparams.cohunatt = 1;
        fixedParams.att = 1;
    end
    fixedParams.fitatt = 1;
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

%% fit HRF
fit = fitModel(data);
if strfind(mode,'useprf')
    fit.useprf = 1;
else
    fit.useprf = 0;
end

function hrfparams = copyhrfparams(fit)
hrfparams.amp1 = fit.params.amp1;
hrfparams.tau1 = fit.params.tau1;
hrfparams.timelag1 = fit.params.timelag1;
hrfparams.amp2 = fit.params.amp2;
hrfparams.tau2 = fit.params.tau2;
hrfparams.timelag2 = fit.params.timelag2;
hrfparams.exponent = fit.params.exponent;
hrfparams.n = fit.params.n;
hrfparams.offset = fit.params.offset;
function roiparams = copyroiparams(fit)
p = fit.params;
f = fields(p);
for fi = 1:length(f)
    cfield = f{fi};
    if any(cellfun(@(x) ~isempty(strfind(cfield,x)),fit.ROIs))
        roiparams.(cfield) = p.(cfield);
    end
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

optimParams = optimset('Algorithm','levenberg-marquardt','MaxIter',inf,'Display','off');
[bestparams, ~, res, ~, ~, ~, curjacob] = lsqnonlin(@hrfResidual,initparams,minparams,maxparams,optimParams,data.tSeries,data.design,data.runtrans,f,fixedParams);

[~,fit] = hrfResidual(bestparams,data.tSeries,data.design,data.runtrans,0,fixedParams);
fit.params = getParams(bestparams,fixedParams);
fit.roiparams = cell(size(fixedParams.ROIs));
for ri = 1:length(fixedParams.ROIs)
    fit.roiparams{ri} = getROIParams(fit.params,fixedParams.ROIs{ri});
end
% fit.roiparams = getROIParams(fit.params,fixedParams.ROIs{ri});
fit.ROIs = fixedParams.ROIs;

function [res, fit] = hrfResidual(params,tSeries,design,runtrans,f,fixedParams)

fit = struct;

    % stimvol basecon lcon rcon basecoh lcoh rcoh timing task
params = getParams(params,fixedParams);

t = 0.25:0.5:49.75;
impulse = gamma(t,params);
if fixedParams.fithrf
    canonImpulse = impulse;
end

res = [];
fit.model = cell(size(fixedParams.ROIs));
fit.tSeries = tSeries;
for ri = 1:length(fixedParams.ROIs)
    roiparams = getROIParams(params,fixedParams.ROIs{ri});
    if fixedParams.fithrf
        impulse = canonImpulse * roiparams.betas;
    end
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
            if fixedParams.fitatt
                % adjust parameters according to current trial
                % attended/unattended
                if cdesign(si,9)==1
                    % attending COHERENCE
                    if fixedParams.att==1
                        roiparams.Rmax = roiparams.Rmax*roiparams.conunatt;
                        roiparams.slope = roiparams.slope*roiparams.cohatt;
                    else
                        warning('Not implemented');
                        keyboard
                        coneff = coneff + roiparams.conunatt;
                        coheff = coheff + roiparams.cohatt;
                    end
                elseif cdesign(si,9)==2
                    % attending CONTRAST
                    if fixedParams.att==1
                        roiparams.Rmax = roiparams.Rmax*roiparams.conatt;
                        roiparams.slope = roiparams.slope*roiparams.cohunatt;
                    else
                        warning('Not implemented');
                        keyboard
                        coneff = coneff + roiparams.conatt;
                        coheff = coheff + roiparams.cohunatt;
                    end
                end
            end
            % okay, for each stimvol, place its effect
            sv = cdesign(si,1);
            if ((cdesign(si,conidx)-cdesign(si,2))+(cdesign(si,cohidx)-cdesign(si,5)))>0
                if fixedParams.fitroi || fixedParams.fitatt
                    coneff = conModel(cdesign(si,conidx)-cdesign(si,2),roiparams);
                    coheff = cohModel(cdesign(si,cohidx)-cdesign(si,5),roiparams);

                    effect = coneff+coheff+params.offset;

                    % adjust timing
                    if cdesign(si,8)>=1
                        idxs = sv:min(runtrans(run,2),sv+cdesign(si,8));
                    else
                        idxs = sv;
                        effect = effect * cdesign(si,8);
                    end
                    if length(effect)>1
                        roimodel(:,idxs) = roimodel(:,idxs)+repmat(effect,1,length(idxs));
                    else
                        roimodel(:,idxs) = roimodel(:,idxs)+repmat(effect,size(roimodel,1),length(idxs));
                    end

                else
                    if cdesign(si,8)>=1
                        idxs = sv:min(runtrans(run,2),sv+cdesign(si,8));
                    else
                        idxs = sv;
                    end
                    roimodel(:,idxs) = 0.5; % why 0.5? Our TR is 500 ms, so 0.5 means that the HRF will be scaled such that it approximates 1% signal change / sec
                end
            end
        end
    end
    cmt = conv(impulse,roimodel);
    fit.model{ri} = cmt(1:size(ctSeries,2));
    res = [res ctSeries-fit.model{ri}];
end

ssres = sum(res.^2);
fit.r2 = 1 - ssres/fixedParams.sstot;

if f>0
    if fixedParams.fithrf
        figure(f)
        clf(f), hold on
        curd = constructD(fit.tSeries{1}/100+1,fixedParams.sv,0.5,25,fixedParams.concatInfo,'none','deconv',0);
        decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
        decon = rmfield(decon,'scm');
        decon = rmfield(decon,'covar');
        plot(decon.ehdr,'o');
        curd = constructD(fit.model{1}/100+1,fixedParams.sv,0.5,25,fixedParams.concatInfo,'none','deconv',0);
        decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
        decon = rmfield(decon,'scm');
        decon = rmfield(decon,'covar');
        plot(decon.ehdr);
        title(sprintf('Amp1: %1.2f Tau1: %1.2f TL1: %1.2f Amp2: %1.2f Tau2: %1.2f TL2: %1.2f ',params.amp1,params.tau1,params.timelag1,params.amp2,params.tau2,params.timelag2));
    else
        figure(f)
        clf(f)
        subplot(211)
        hold on
        plot(ctSeries(1:1000),'b');
        plot(fit.model{end}(1:1000),'r');
        title(sprintf('R^2: %0.2f',fit.r2));
        subplot(212)
        plot(t,impulse);
    end
end

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

function [initparams, minparams, maxparams] = initParams()
%%
global fixedParams params

%% Deal with HRF params
fixedParams.strs = fields(params.hrfparams)';
rfields = fields(params.roiparams);

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