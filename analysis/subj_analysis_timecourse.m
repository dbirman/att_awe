function subj_analysis_timecourse( subj, incl )

ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','LO1','LO2','MT'};


%% Joint Analysis for s0025

%% Load Files
fdir = '~/data/cohcon_localizer'; 
files = dir(fullfile(fdir,sprintf('%s*.mat',subj)));


analyses = {'CohxCon','Timing'};
%% Load each file
datas = {};

for fi = 1:length(files)
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
        rts{end+1} = datas{fi}.pre.tSeries;
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


%% Plot
for ri = 1:length(bestfits)
    try
    bestfit = bestfits{ri};
    
    clist = brewermap(20,'PuOr');
    res = [0.25 0.5 1 2 4];
%     if ~isdir(fullfile(pwd,'Figures/linear')), mkdir(fullfile(pwd,'Figures/linear')); end

    figure, hold on
    
    count = 1;
    for i = incl
        subplot(length(datas),1,i), hold on
        for j = 1:size(datas{i}.deconvo{ri}.ehdr,1)
            % plot original
            plot(time,datas{i}.deconvo{ri}.ehdr(j,:),'o','MarkerFaceColor',clist(j,:),'MarkerEdgeColor',[1 1 1]);
            h = errbar(time,datas{i}.deconvo{ri}.ehdr(j,:),datas{i}.deconvo{ri}.ehdrste(j,:),'Color',clist(j,:));
        end
        
        for k = 1:size(datas{i}.deconvo{ri}.ehdr,1)
            % plot the fit
            plot(time2,bestfit.out(count,:),'-','Color',clist(k,:));
            count = count + 1;
        end
        title(analyses{i})
    end
    axis([0 30 -0.5 1.5]);
    xlabel('Time (s)');
    ylabel('BOLD Amplitude (%s Signal Change)');
    title(sprintf('25%% Coherence, R^2: %0.2f',bestfit.r2));
        drawPublishAxis;

    %%
    figure, hold on
    plot(bestfit.full.fcon,bestfit.full.fconr,'r');
    plot(bestfit.full.fcoh,bestfit.full.fcohr,'b');
    title(ROIs{ri});
    %%
    if ~isdir(fullfile(pwd,'Figures')), mkdir(fullfile(pwd,'Figures')); end
    if ~isdir(fullfile(pwd,'Figures/linear')), mkdir(fullfile(pwd,'Figures/linear')); end
    fname = fullfile(pwd,'Figures/linear',sprintf('%s_%s_lintiming.pdf',subj,ROIs{ri}));
    print(fname,'-dpdf')
    catch
        
    if ~isdir(fullfile(pwd,'Figures')), mkdir(fullfile(pwd,'Figures')); end
    if ~isdir(fullfile(pwd,'Figures/linear')), mkdir(fullfile(pwd,'Figures/linear')); end
    fname = fullfile(pwd,'Figures/linear',sprintf('%s_%s_lintiming.pdf',subj,ROIs{ri}));
    print(fname,'-dpdf')
    end
end

%% Compute Effect
effect = struct;
effect.con = [];
effect.coh = [];
effect.ROIs = ROIs;

for ri = 1:length(bestfits)
    bestfit = bestfits{ri};
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