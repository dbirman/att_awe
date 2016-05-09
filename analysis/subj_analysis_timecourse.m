function subj_analysis_timecourse( subj, incl )
 
ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','LO1','LO2','MT'};

%% Load Files
fdir = '~/data/cohcon_localizer'; 
files = dir(fullfile(fdir,sprintf('%s2016*.mat',subj)));


analyses = {'CohxCon','Timing'};
%% Load each file
datas = {};

for fi = 1
    load(fullfile(fdir,files(fi).name));
    datas{fi} = data;
end

%% Move everything into these convenient holders

timeseries = {};
runtrans = {};
stimvol = {};
stimnames = {};
basecon = {};
basecoh = {};
con = {};
coh = {};
timing = {};

for ri = 1:length(ROIs)
    rts = {};
    rrt = {};
    rsv = {};
    rsn = {};
    rbcn = {};
    rbch = {};
    rcn = {};
    rch = {};
    rt = {};
    for fi = 1:length(datas)
        rts{end+1} = datas{fi}.pre.tSeries{ri};
        rrt{end+1} = datas{fi}.pre.concatInfo.runTransition;
        rsv{end+1} = datas{fi}.pre.stimvol;
        rsn{end+1} = datas{fi}.pre.stimNames;
        rbcn{end+1} = datas{fi}.basecon;
        rbch{end+1} = datas{fi}.basecoh;
        rcn{end+1} = datas{fi}.con;
        rch{end+1} = datas{fi}.coh;
        rt{end+1} = datas{fi}.timing;
    end
    timeseries{end+1} = rts;
    runtrans{end+1} = rrt;
    stimvol{end+1} = rsv;
    stimnames{end+1} = rsn;
    basecon{end+1} = rbcn;
    basecoh{end+1} = rbch;
    con{end+1} = rcn;
    coh{end+1} = rch;
    timing{end+1} = rt;
end

%% Fit Timecourse model
fits = cell(size(ROIs));
for ri = 1:length(ROIs)
    fits{ri} = fitCCTimecourseVoxelModel(timeseries{ri},stimvol{ri},runtrans{ri},basecon{ri},basecoh{ri},stimnames{ri},timing{ri});
end

%% Check impulses
clist = brewermap(10,'YlOrBr');
figure, hold on
for ri = 1:length(fits)
    plot(fits{ri}.impulse,'Color',clist(ri,:));
end
legend(ROIs)

%% Con/Coh model plots

% we want a plot that shows the change in the contrast/coherence functions
% across brain regions, for now just show v1-v2-v3-v3a-mt
roiIndexes = [1 2 3 5 10];
clist = brewermap(5,'Oranges');
figure
subplot(211), hold on
for ri = roiIndexes
    plot(fits{ri}.full.fcon,fits{ri}.full.fconr,'Color',clist(find(roiIndexes==ri,1),:));
end
legend(ROIs(roiIndexes),'Location','NorthWest');
title('Contrast Sensitivity by ROI')
xlabel('Contrast (%)');
ylabel('Effect (% signal change / s)');
drawPublishAxis
subplot(212), hold on
clist = brewermap(5,'Purples');
for ri = roiIndexes
    plot(fits{ri}.full.fcoh,fits{ri}.full.fcohr,'Color',clist(find(roiIndexes==ri,1),:));
end
legend(ROIs(roiIndexes),'Location','NorthWest');
title('Coherence Sensitivity by ROI')
xlabel('Motion Coherence (%)');
ylabel('Effect (% signal change / s)');
drawPublishAxis

if ~isdir(fullfile(pwd,'Figures')), mkdir(fullfile(pwd,'Figures')); end
if ~isdir(fullfile(pwd,'Figures/sensitivity')), mkdir(fullfile(pwd,'Figures/sensitivity')); end
fname = fullfile(pwd,'Figures/sensitivity',sprintf('%s_sensitivity.pdf',subj));
print(fname,'-dpdf')

%% Deconvolve conditions

for ri = 10
    decon = cell(length(fits),length(datas));
    model = cell(length(fits),length(datas));
    for di = 1:length(datas)
        curd = constructD(fits{ri}.orig.timeseries{di}-1,stimvol{ri}{di},0.5,15,datas{di}.pre.concatInfo,'none','deconv',0);

        decon{ri,di} = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
        decon{ri,di} = rmfield(decon{ri,di},'scm');
        decon{ri,di} = rmfield(decon{ri,di},'covar');
        
        modd = constructD(fits{ri}.out{di},stimvol{ri}{di},0.5,15,datas{di}.pre.concatInfo,'none','deconv',0);

        model{ri,di} = getr2timecourse(modd.timecourse,modd.nhdr,modd.hdrlenTR,modd.scm,modd.framePeriod,modd.verbose);
        model{ri,di} = rmfield(model{ri,di},'scm');
        model{ri,di} = rmfield(model{ri,di},'covar');
    end
end

%% Plot Fit Comparisons
for ri = 10
    figure
    for di = 1:length(datas)
        subplot(length(datas),1,di), hold on
        
        [con,coh,time] = parseNames(stimnames{ri}{di},'contrast=','coherence=','timing=',' and ');
        d = decon{ri,di};
        m = model{ri,di};
        
        if isempty(time)
            %coh x con
            % how best to display this? for now just do the base conditions
            % where either con/coh stayed constant
            const_con = find(con==0.25);
            const_coh = find(coh==0);
            clist = brewermap(11,'PuOr');
            for i = 1:5 % coherences
                plot(d.time,d.ehdr(const_con(i),:),'o','MarkerFaceColor',clist(11-i,:),'MarkerEdgeColor',[1 1 1]);
                plot(m.time,m.ehdr(const_con(i),:),'-','Color',clist(11-i,:));
            end
            for i = 1:4 % contrasts
                plot(d.time,d.ehdr(const_coh(i),:),'o','MarkerFaceColor',clist(i,:),'MarkerEdgeColor',[1 1 1]);
                plot(m.time,m.ehdr(const_coh(i),:),'-','Color',clist(i,:));
            end
        else
            %timing
            clist = brewermap(5,'Purples');
            for i = 1:5
                plot(d.time,d.ehdr(i,:),'o','MarkerFaceColor',clist(i,:),'MarkerEdgeColor',[1 1 1]);
                plot(m.time,m.ehdr(i,:),'-','Color',clist(i,:));
            end
            for i = 6:10
                plot(d.time,d.ehdr(i,:),'o','MarkerFaceColor',clist(i-5,:),'MarkerEdgeColor',[1 1 1]);
                plot(m.time,m.ehdr(i,:),'-','Color',clist(i-5,:));
            end
        end
    end
end

%% Compute Effect
effect = struct;
effect.con = [];
effect.coh = [];
effect.ROIs = ROIs;

for ri = 1:length(fits)
    bestfit = fits{ri};
    effect.con = [effect.con bestfit.full.fconr(end)];
    effect.coh = [effect.coh bestfit.full.fcohr(end)];
end

% save for cross subjecte ffects
fname = fullfile(pwd,sprintf('%s_effect.mat',subj));
save(fname,'effect')

%% Compute Effect Plot
if ~isdir(fullfile(pwd,'Figures/effect')),mkdir(fullfile(pwd,'Figures/effect')); end
fname = fullfile(pwd,'Figures/effect',sprintf('%s_effect.pdf',subj));

figure, hold on
cmap = brewermap(2,'PuOr');
spacer = 1;
for ri = 1:length(ROIs)
    if (effect.con(ri)<=3) && (effect.coh(ri)<=4)
        bar(spacer*ri-0.2,effect.con(ri),0.25,'FaceColor',cmap(1,:));
        bar(spacer*ri+0.2,effect.coh(ri),0.25,'FaceColor',cmap(2,:));
    end
%     errbar(spacer*ri-0.2,coneff(ri,1),coneffint(ri,1),'Color',[0 0 0]);
%     errbar(spacer*ri+0.2,coheff(ri,1),coheffint(ri,1),'Color',[0 0 0]);
end
set(gca,'XTick',(1:length(ROIs))*spacer,'XTickLabel',ROIs);
ylabel('Sensitivity (Regression Slope)');
legend({'Contrast','Coherence'});
drawPublishAxis
print(fname,'-dpdf')