%% Load all data and average across shared ROIs
pfxs = {'l','r'};
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

allROIs = {};
for ri = 1:length(rois)
    for pi = 1:2
        allROIs{end+1} = sprintf('%s%s',pfxs{pi},rois{ri});
    end
end

% load each subject's relevant data and collapse into a giant monster
% array
allCon = zeros(length(nSIDs),16,16,81);
allCoh = zeros(length(nSIDs),16,16,81);
for ni = 1:length(nSIDs)
    fname = fullfile(datafolder,sprintf('s%04.0f_deconEffects_att.mat',nSIDs(ni)));
    load(fname);
    
    for ri = 1:length(allROIs)
        allCon(ni,ri,:,:) = decondata.(allROIs{ri}).cc.conresp;
        allCoh(ni,ri,:,:) = decondata.(allROIs{ri}).cc.cohresp;
    end
end
    
% average across hemispheres
allCon_ = zeros(length(nSIDs),8,16,81);
allCoh_ = zeros(length(nSIDs),8,16,81); 
for ri = 1:length(rois)
    roinums = cellfun(@(x) ~isempty(strfind(x,rois{ri})),allROIs,'UniformOutput',false);
    roinums = find([roinums{:}]);
    
    if length(roinums)>2 % V3 being mis-interpreted as V3a/V3b
        roinums = roinums(1:2);
    end
    
    allCon_(:,ri,:,:) = squeeze(mean(allCon(:,roinums,:,:),2));
    allCoh_(:,ri,:,:) = squeeze(mean(allCoh(:,roinums,:,:),2));
end

% average across subjects

mCon = squeeze(bootci(1000,@nanmean,allCon_));
mCoh = squeeze(bootci(1000,@nanmean,allCoh_));
mCon_ = squeeze(mean(mCon));
mCoh_ = squeeze(mean(mCoh));

%% Plot each set of conditions
conidx = decondata.lV1.cc.conidxs;
cohidx = decondata.lV1.cc.cohidxs;
taskidx = decondata.lV1.cc.taskidx;
deltaidx = decondata.lV1.cc.deltaidx;
conStim = decondata.lV1.cc.conStim;
cohStim = decondata.lV1.cc.cohStim;
contrasts = [0.325 0.4 0.55 0.85]-.25;
coherences = [0.15 0.3 0.45 0.6];

%% Prep save
avgdecon.allCon = allCon_;
avgdecon.allCoh = allCoh_;
avgdecon.mCon = mCon; % includes confidence intervals
avgdecon.mCoh = mCoh;
avgdecon.conidx = conidx;
avgdecon.cohidx = cohidx;
avgdecon.taskidx = taskidx;
avgdecon.deltaidx = deltaidx;
avgdecon.conStim = conStim;
avgdecon.cohStim = cohStim;
avgdecon.contrasts = contrasts;
avgdecon.coherences = coherences;

%% Remove delta==1
% conidx = conidx(logical(~deltaidx));
% taskidx = taskidx(logical(~deltaidx));
% constim = conStim(logical(~deltaidx));
% cohstim = cohStim(logical(~deltaidx));


%% CONTRAST
for ri = 1:8
    h = figure; hold on

    taskdisp = {'-o','--o'};
    concolor = brewermap(11,'PuOr');
    concolor = flipud(concolor(1:4,:));

    resp = squeeze(mCon_(ri,:,:));
%     resp = resp(logical(~deltaidx),:);

%     stimNames = stimNames(logical(~deltaidx));

    for si = 1:size(resp,1)
        subplot(1,2,taskidx(si)); hold on
        plot(0.25:.5:40.25,resp(si,:),'o','MarkerSize',10,'MarkerFaceColor',concolor(find(conidx(si)==contrasts,1),:),'MarkerEdgeColor',[1 1 1]);
    end

    tasks = {'Motion','Contrast'};
    for i = 1:2
        subplot(1,2,i);
        title(sprintf('%s: Task %s',rois{ri},tasks{i}));
        a = axis;
        axis([0 15 -1 1.5]);
        set(gca,'YTick',[-0.5 0 0.5 1]);

%         if i==2, legend(constim); end
        drawPublishAxis
    end
    savepdf(h,fullfile(datafolder,'avg_deconatt',sprintf('contrast_%s.pdf',ROIs{ri})));
end
%% MOTION 
for ri = 1:8
    h = figure; hold on

    taskdisp = {'-o','--o'};
    cohcolor = brewermap(11,'PuOr');
    cohcolor = cohcolor(8:11,:);


    resp = squeeze(mCoh_(ri,:,:));
    resp = resp(logical(~deltaidx),:);

%     stimNames = stimNames(logical(~deltaidx));

    for si = 1:size(resp,1)
        subplot(1,2,taskidx(si)); hold on
        plot(0.25:.5:40.25,resp(si,:),'o','MarkerSize',10,'MarkerFaceColor',cohcolor(find(cohidx(si)==coherences,1),:),'MarkerEdgeColor',[1 1 1]);
    end

    tasks = {'Motion','Contrast'};
    for i = 1:2
        subplot(1,2,i);
        title(sprintf('%s: Task %s',rois{ri},tasks{i}));
        a = axis;
        axis([0 15 -0.5 1.25]);
        set(gca,'YTick',[-0.5 0 0.5 1]);

%         if i==2, legend(stimNames); end
        drawPublishAxis
    end
    savepdf(h,fullfile(datafolder,'avg_deconatt',sprintf('motion_%s.pdf',ROIs{ri})));
end

%% Save

fname = fullfile(datafolder,'avg_deconatt','avg_deconEffects_att.mat');
save(fname,'avgdecon');