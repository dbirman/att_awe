function deconCohCon( data, fit, roinum )
%DECONCOHCON Deconvolve the cohxcon experiment

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
            if fit.useprf
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

%%
% Split the design into TIMING, CONTRAST_x_COHERENCE, and TASK
% This assumes that the design is already concatenated (i.e. can't deal
% with {datas} structure

task_idxs = logical(~(data.design(:,9)==0));
cohxcon_idxs = logical(logical(data.design(:,8)==5).*logical(~task_idxs));
timing_idxs = logical((data.design(:,8)~=5).*~task_idxs);

taskdesign = data.design(task_idxs,:);
cohxcondesign = data.design(cohxcon_idxs,:);
timingdesign = data.design(timing_idxs,:);

designs = {cohxcondesign,timingdesign,taskdesign};

%% reduce space of cohxcon values (for task)

%% deconvolve for each timecourse

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
concatInfo.runTransition = data.runtrans;
curd = constructD(data.tSeries{roinum},timing_sv,0.5,30,concatInfo,'none','deconv',0);
decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
decon = rmfield(decon,'scm');
decon = rmfield(decon,'covar');
curd = constructD(fit.model{roinum}/100+1,timing_sv,0.5,30,concatInfo,'none','deconv',0);
model = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
model = rmfield(model,'scm');
model = rmfield(model,'covar');

%% Time plot
clist = brewermap(30,'PuOr');
figure
% 5 plots
subplot(4,2,1:4), hold on
title('Contrast +0%');
idxs = find(conidxs==0.25);
for i = idxs
    plot(decon.time,decon.ehdr(i,:),'o','MarkerFaceColor',clist(i,:),'MarkerEdgeColor',[1 1 1]);
    errbar(decon.time,decon.ehdr(i,:),decon.ehdrste(i,:),'Color',clist(i,:));
    plot(model.time,model.ehdr(i,:),'Color',clist(i,:));
end

subplot(4,2,5), hold on
title('Contrast +25% Coherence +0%');
idxs = find((conidxs==0.5).*(cohidxs==0.25));
for i = idxs
    plot(decon.time,decon.ehdr(i,:),'o','MarkerFaceColor',clist(i,:),'MarkerEdgeColor',[1 1 1]);
    errbar(decon.time,decon.ehdr(i,:),decon.ehdrste(i,:),'Color',clist(i,:));
    plot(model.time,model.ehdr(i,:),'Color',clist(i,:));
end

subplot(4,2,6), hold on
title('Contrast +25% Coherence +100%');
idxs = find((conidxs==0.5).*(cohidxs==1));
for i = idxs
    plot(decon.time,decon.ehdr(i,:),'o','MarkerFaceColor',clist(i,:),'MarkerEdgeColor',[1 1 1]);
    errbar(decon.time,decon.ehdr(i,:),decon.ehdrste(i,:),'Color',clist(i,:));
    plot(model.time,model.ehdr(i,:),'Color',clist(i,:));
end
subplot(4,2,7), hold on
title('Contrast +75% Coherence +0%');
idxs = find((conidxs==1).*(cohidxs==0.25));
for i = idxs
    plot(decon.time,decon.ehdr(i,:),'o','MarkerFaceColor',clist(i,:),'MarkerEdgeColor',[1 1 1]);
    errbar(decon.time,decon.ehdr(i,:),decon.ehdrste(i,:),'Color',clist(i,:));
    plot(model.time,model.ehdr(i,:),'Color',clist(i,:));
end
subplot(4,2,8), hold on
title('Contrast +75% Coherence +100%');
idxs = find((conidxs==1).*(cohidxs==1));
for i = idxs
    plot(decon.time,decon.ehdr(i,:),'o','MarkerFaceColor',clist(i,:),'MarkerEdgeColor',[1 1 1]);
    errbar(decon.time,decon.ehdr(i,:),decon.ehdrste(i,:),'Color',clist(i,:));
    plot(model.time,model.ehdr(i,:),'Color',clist(i,:));
end
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
concatInfo.runTransition = data.runtrans;
curd = constructD(data.tSeries{roinum},cohxcon_sv,0.5,30,concatInfo,'none','deconv',0);
decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
decon = rmfield(decon,'scm');
decon = rmfield(decon,'covar');
curd = constructD(fit.model{roinum}/100+1,cohxcon_sv,0.5,30,concatInfo,'none','deconv',0);
model = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
model = rmfield(model,'scm');
model = rmfield(model,'covar');

clist = brewermap(10,'PuOr');
figure
% plot #1: contrast against 0% coherence (no change)
lconidx = find(cohidx==0);
subplot(2,1,1), hold on
for i = 1:4
    plot(decon.time,decon.ehdr(lconidx(i),:),'o','MarkerFaceColor',clist(i,:),'MarkerEdgeColor',[1 1 1]);
    errbar(decon.time,decon.ehdr(lconidx(i),:),decon.ehdrste(lconidx(i),:),'Color',clist(i,:));
end
for i = 1:4
    plot(model.time,model.ehdr(lconidx(i),:),'Color',clist(i,:));
end
legend({'Contrast=+0%','Contrast=+25%','Contrast=+50%','Contrast=+75%'});
% plot #2: coherence against 25% contrast (no change)
lcohidx = find(conidx==0.25);
subplot(2,1,2), hold on
for i=1:5
    plot(decon.time,decon.ehdr(lcohidx(i),:),'o','MarkerFaceColor',clist(i+5,:),'MarkerEdgeColor',[1 1 1]);
    errbar(decon.time,decon.ehdr(lcohidx(i),:),decon.ehdrste(lcohidx(i),:),'Color',clist(i+5,:));
end
for i = 1:5
    plot(model.time,model.ehdr(lcohidx(i),:),'Color',clist(i+5,:));
end
legend({'Coherence+0%','Coherence+25%','Coherence+50%','Coherence+75%','Coherence+100%'});
