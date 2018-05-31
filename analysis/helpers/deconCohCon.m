function [rr2, sd] = deconCohCon( data, roiname ,subj)
%DECONCOHCON Deconvolve the cohxcon experiment

% 6/17 update: adding cross-validation. 

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
    data.rtSeriesAll = cell(1,length(data_old{1}.ROIs));
    data.rtSeries25 = cell(1,length(data_old{1}.ROIs));
    data.rtSeries2 = cell(1,length(data_old{1}.ROIs));
    length_sofar = 0;
    for ri = 1:length(data.tSeries), data.tSeries{ri} = []; end
    data.design = [];
    data.runtrans = [];
    for si = 1:length(data_old)
        for ri = 1:length(data_old{si}.tSeries)
            data.tSeries{ri} = [data.tSeries{ri} data_old{si}.tSeries{ri}];
            data.rtSeriesAll{ri} = [data.rtSeriesAll{ri} data_old{si}.rtSeriesAll{ri}];
            data.rtSeries25{ri} = [data.rtSeries25{ri} data_old{si}.rtSeries25{ri}];
            data.rtSeries2{ri} = [data.rtSeries2{ri} data_old{si}.rtSeries2{ri}];
        end
        % tweak the SV by adding 
        cdes = data_old{si}.design;
        cdes(:,1) = cdes(:,1) + length_sofar;
        data.design = [data.design ; cdes];
        data.runtrans = [data.runtrans ; data_old{si}.runtrans+length_sofar];
        length_sofar = length_sofar + length(data_old{si}.tSeries{1});
    end
end

%%
% if ~isempty(strfind(fit.mode,'droptiming'))
%     disp('%% All Timing Responses Dropped %%');
%     idxs = data.design(:,8)==5;
%     data.design = data.design(idxs,:);
% end

%% Concatenate multiple ROIs when requested
roinums = cellfun(@(x) ~isempty(strfind(x,roiname)),data.ROIs,'UniformOutput',false);
roinums = find([roinums{:}]);

if length(roinums)>2
    disp('Found more than two. Probably V3 found V3a/V3b. Taking first two.');
    roinums = roinums(1:2);
end

if length(roinums)==2
    disp('Using two ROIs concatenated');
    tSeries = [data.tSeries{roinums(1)} data.tSeries{roinums(2)}];
    rtSeriesAll = [data.rtSeriesAll{roinums(1)} data.rtSeriesAll{roinums(2)}];
    rtSeries25 = [data.rtSeries25{roinums(1)} data.rtSeries25{roinums(2)}];
    rtSeries2 = [data.rtSeries2{roinums(1)} data.rtSeries2{roinums(2)}];
    % double the design
    design2 = data.design;
    design2(:,1) = design2(:,1)+length(data.tSeries{1}); % add to stimvols
    design = [data.design ; design2];
    % run transitions
    runtrans2 = data.runtrans+length(data.tSeries{1});
    runtrans = [data.runtrans ; runtrans2];
else
    disp('failure');
    keyboard
    roinum = roinums(1);
    tSeries = data.tSeries{roinum};
    design = data.design;
    runtrans = data.runtrans;
end

%%
% Split the design into TIMING, CONTRAST_x_COHERENCE, and TASK
% This assumes that the design is already concatenated (i.e. can't deal
% with {datas} structure

% [csv basecon lCon rCon basecoh lCoh rCoh tim task];
% for task, remove everything where no task was given (fixation task)
task_idxs = logical(~(design(:,9)==0));
% for cohcon, remove anything where the time differs from 2.5 s
cohxcon_idxs = logical(logical(design(:,8)==5).*logical(~task_idxs));
% for timing, take what's left
timing_idxs = logical((design(:,8)~=5).*~task_idxs);

taskdesign = design(task_idxs,:);
cohxcondesign = design(cohxcon_idxs,:);
timingdesign = design(timing_idxs,:);

designs = {cohxcondesign,timingdesign,taskdesign};

%% reduce space of cohxcon values (for task)


%% CohxCon
contrast = cohxcondesign(:,3); ucon = unique(contrast);
coherence = cohxcondesign(:,6); ucoh = unique(coherence);
sv = cohxcondesign(:,1);
cohxcon_sv = {};
stimNames = {};
conidx = [];
cohidx = [];
for coni = 1:length(ucon)
    for cohi = 1:length(ucoh)
        cohxcon_sv{end+1} = sv(logical((contrast==ucon(coni)).*(coherence==ucoh(cohi))));
        stimNames{end+1} = sprintf('Contrast=%0.2f Coherence=%0.2f',ucon(coni),ucoh(cohi));
        conidx(end+1) = ucon(coni);
        cohidx(end+1) = ucoh(cohi);
    end
end
concatInfo.runTransition = runtrans;

%% Load previous cohxcon data

fname = fullfile(datafolder,sprintf('%s_decon.mat',subj));
if exist(fname,'file')==2, load(fname); end

%% Deconvolve
tsOpts = {'tSeries','rtSeriesAll','rtSeries2','rtSeries25'};
saveNames = {'resp_prf','resp_all','resp_2','resp_25'};

for ti = 1:length(tsOpts)
    cts = eval(tsOpts{ti});
    curd = constructD(cts,cohxcon_sv,0.5,40,concatInfo,'none','deconv',0);
    decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
    decon = rmfield(decon,'scm');
    decon = rmfield(decon,'covar');
    temp = std(decon.resid{1,1});
    sd(ti,1) = temp(1);
%     if ti==4
        rr2(ti,1) = decon.r2;
%     end
    decondata.(roiname).cc.(saveNames{ti}) = decon.ehdr;
end
% curd = constructD(mtSeries/100+1,cohxcon_sv,0.5,40,concatInfo,'none','deconv',0);
% model = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
% model = rmfield(model,'scm');
% model = rmfield(model,'covar');

%% Save cohxcon data
disp('not saving');
% decondata.(roiname).cc.conidxs = conidx;
% decondata.(roiname).cc.cohidxs = cohidx;
% save(fname,'decondata');

%% Return if no timing required
% if ~isempty(strfind(mode,'droptiming'))
%     return
% end
% return
%% Timing
%  - Ignore l vs. right since identical
contrast = timingdesign(:,3); ucon = unique(contrast);
coherence = timingdesign(:,6); ucoh = unique(coherence);
timing = timingdesign(:,8); ut = unique(timing);
sv = timingdesign(:,1);
timing_sv = {};
conidxs = [];
cohidxs = [];
timidxs = [];
for coni = 1:length(ucon)
    for cohi = 1:length(ucoh)
        for ti = 1:length(ut)
            timing_sv{end+1} = sv(logical((contrast==ucon(coni)).*(coherence==ucoh(cohi)).*(timing==ut(ti))));
            conidxs(end+1) = ucon(coni);
            cohidxs(end+1) = ucoh(cohi);
            timidxs(end+1) = ut(ti);
        end
    end
end

%% Load previous data

fname = fullfile(datafolder,sprintf('%s_decon.mat',subj));
if exist(fname,'file')==2, load(fname); end

%% Deconvolve
tsOpts = {'tSeries','rtSeriesAll','rtSeries2','rtSeries25'};
saveNames = {'resp_prf','resp_all','resp_2','resp_25'};

for ti = 1:length(tsOpts)
    cts = eval(tsOpts{ti});
    curd = constructD(cts,timing_sv,0.5,40,concatInfo,'none','deconv',0);
    decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
    decon = rmfield(decon,'scm');
    decon = rmfield(decon,'covar');
    temp = std(decon.resid{1,1});
    sd(ti,2) = temp(1);
%     if ti==4
        rr2(ti,2) = decon.r2;
%     end
    decondata.(roiname).time.(saveNames{ti}) = decon.ehdr;
end

%% Save timing data
disp('not saving');
% decondata.(roiname).time.conidxs = conidxs;
% decondata.(roiname).time.cohidxs = cohidxs;
% decondata.(roiname).time.timidxs = timidxs;
% save(fname,'decondata');
