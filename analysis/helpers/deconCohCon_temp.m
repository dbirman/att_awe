function deconCohCon( data, roiname ,subj)
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

%% Re-discover runs, because you're an idiot and you didn't build CV in from the start
% Using designs, find which runs are cohcon and which are timing. Runtrans
% has the actual timepoints corresponding to the runs
cohcon_runs = [];
timing_runs = [];

for ri = 1:size(runtrans,1)
    cc_sv = cohxcondesign(:,1);
    cc_sv = cc_sv(cc_sv>=runtrans(ri,1));
    cc_sv = cc_sv(cc_sv<=runtrans(ri,2));
    if ~isempty(cc_sv)
        cohcon_runs(end+1) = ri;
    else
        timing_runs(end+1) = ri;
    end
end

%% CohxCon
tsOpts = {'tSeries','rtSeriesAll','rtSeries2','rtSeries25'};
contrast = cohxcondesign(:,3); ucon = unique(contrast);
coherence = cohxcondesign(:,6); ucoh = unique(coherence);

cc = struct;
for ti = 1:4
    disp(sprintf('Running time series %s ',tsOpts{ti}));
    ctSeries = tsOpts{ti};
    curts = eval(ctSeries);
    disppercent(-1/length(cohcon_runs));
    for ci = 1:length(cohcon_runs)
        fold = sprintf('fold%i',ci);
        disp(fold);
        cc.(fold) = struct;
        cc.(fold).(ctSeries) = struct;
        % remove stimvols that are part of the test set for this run
        % compute the test set
        svgroup = logical((cohxcondesign(:,1)>=runtrans(cohcon_runs(ci),1)) .* (cohxcondesign(:,1)<=runtrans(cohcon_runs(ci),2)));
        testdesign = cohxcondesign(svgroup,:);
        traindesign = cohxcondesign(~svgroup,:);

        sv_train = traindesign(:,1);
        sv_test = testdesign(:,1);

        sv_train_ = {};
        sv_test_ = {};

        stimNames = {};
        conidx = [];
        cohidx = [];
        for coni = 1:length(ucon)
            for cohi = 1:length(ucoh)
                sv_train_{end+1} = sv_train(logical((traindesign(:,3)==ucon(coni)).*((traindesign(:,6)==ucoh(cohi)))));
                sv_test_{end+1} = sv_test(logical(((testdesign(:,3)==ucon(coni)).*((testdesign(:,6)==ucoh(cohi))))));
                stimNames{end+1} = sprintf('Contrast=%0.2f Coherence=%0.2f',ucon(coni),ucoh(cohi));
                conidx(end+1) = ucon(coni);
                cohidx(end+1) = ucoh(cohi);
            end
        end
        concatInfo.runTransition = runtrans;
        traind = constructD(curts,sv_train_,0.5,40,concatInfo,'none','deconv',0);
        decontrain = getr2timecourse(traind.timecourse,traind.nhdr,traind.hdrlenTR,traind.scm,traind.framePeriod,traind.verbose);
        decontrain = rmfield(decontrain,'scm');
        decontrain = rmfield(decontrain,'covar');
        testd = constructD(curts,sv_test_,0.5,40,concatInfo,'none','deconv',0);
        decontest = getr2timecourse(testd.timecourse,testd.nhdr,testd.hdrlenTR,testd.scm,testd.framePeriod,testd.verbose);
        decontest = rmfield(decontest,'scm');
        decontest = rmfield(decontest,'covar');
        cc.(fold).conidxs = conidx;
        cc.(fold).cohidxs = cohidx;
        cc.(fold).(ctSeries).trainresp = decontrain.ehdr;
        cc.(fold).(ctSeries).testresp = decontest.ehdr;
        disppercent(ci/length(cohcon_runs));
    end
    disppercent(inf);
end

%% save cohxcon data
fname = fullfile(datafolder,sprintf('%s_decon.mat',subj));
if exist(fname,'file')==2, load(fname); end
decondata.(roiname).cc = cc;
save(fname,'decondata');
%%
clist = brewermap(10,'PuOr');
f= figure;
% plot #1: contrast against 0% coherence (no change)
lconidx = find(cohidx==0);
convalues = [0.25 0.5 0.75 1];
colmap = [4 3 2 1];
subplot(2,1,1), hold on
for i = 1:4
    ci = colmap(find(convalues==conidx(lconidx(i)),1));
    plot(decon.time,decon.ehdr(lconidx(i),:),'o','MarkerSize',10,'MarkerFaceColor',clist(ci,:),'MarkerEdgeColor',[1 1 1]);
    errbar(decon.time,decon.ehdr(lconidx(i),:),decon.ehdrste(lconidx(i),:),'Color',clist(ci,:));
end
for i = 1:4
    ci = colmap(find(convalues==conidx(lconidx(i)),1));
    plot(model.time,model.ehdr(lconidx(i),:),'Color',clist(ci,:));
end
% legend({'Contrast=+0%','Contrast=+25%','Contrast=+50%','Contrast=+75%'});
% plot #2: coherence against 25% contrast (no change)
lcohidx = find(conidx==0.25);
a = axis;
axis([0 14 a(3) a(4)]);
drawPublishAxis
subplot(2,1,2), hold on
for i=1:5
    plot(decon.time,decon.ehdr(lcohidx(i),:),'o','MarkerSize',10,'MarkerFaceColor',clist(i+5,:),'MarkerEdgeColor',[1 1 1]);
    errbar(decon.time,decon.ehdr(lcohidx(i),:),decon.ehdrste(lcohidx(i),:),'Color',clist(i+5,:));
end
for i = 1:5
    plot(model.time,model.ehdr(lcohidx(i),:),'Color',clist(i+5,:));
end
% legend({'Coherence+0%','Coherence+25%','Coherence+50%','Coherence+75%','Coherence+100%'});

a = axis;
axis([0 14 a(3) a(4)]);
drawPublishAxis
%% Print Figure #2

fname = fullfile(datafolder,sprintf('%s_%s_cohconplot.pdf',subj,roiname));
savepdf(f,fname);

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
concatInfo.runTransition = runtrans;
curd = constructD(tSeries,timing_sv,0.5,40,concatInfo,'none','deconv',0);
decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
decon = rmfield(decon,'scm');
decon = rmfield(decon,'covar');
curd = constructD(mtSeries/100+1,timing_sv,0.5,40,concatInfo,'none','deconv',0);
model = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
model = rmfield(model,'scm');
model = rmfield(model,'covar');

%% Save timing data
fname = fullfile(datafolder,sprintf('%s_decon.mat',subj));
if exist(fname,'file')==2, load(fname); end
decondata.(roiname).time.conidxs = conidxs;
decondata.(roiname).time.cohidxs = cohidxs;
decondata.(roiname).time.timidxs = timidxs;
decondata.(roiname).time.resp = decon.ehdr;
decondata.(roiname).time.mresp = model.ehdr;
save(fname,'decondata');
%% Time plot
clist = brewermap(5,'Greys');
f = figure; hold on
% first plot will be contrast timing, when contrast goes up to 50% or 100%,
% we'll draw each of these with increasing contrast colors
conopts = [0.50 0.50 1 1];
cohopts = [0.25 1 0.25 1];
flip = [0.5 1 2 4 8];
colpos = [1 2 3 4 5];
for sub = 1:4
    subplot(2,2,sub); hold on
    title(sprintf('Con: %i%% Coh: %i%%',conopts(sub)*100,cohopts(sub)*100));
    idxs = find(logical(conidxs==conopts(sub)).*logical(cohidxs==cohopts(sub)));
    for i = idxs
        ci = colpos(find(flip==timidxs(i),1));
        plot(decon.time,decon.ehdr(i,:),'o','MarkerSize',10,'MarkerFaceColor',clist(ci,:),'MarkerEdgeColor',[1 1 1]);
        errbar(decon.time,decon.ehdr(i,:),decon.ehdrste(i,:),'Color',clist(ci,:));
        plot(model.time,model.ehdr(i,:),'Color',clist(ci,:));
    end
    xlabel('Time (s)');
    ylabel('Response (%signal/s)');
    axis([0 15 -2 5]);
    drawPublishAxis
end

%% Print Figure #1
title(sprintf('%s: %s',subj,roiname));
savepdf(f,fullfile(datafolder,sprintf('%s_%s_timeplot.pdf',subj,roiname)));
