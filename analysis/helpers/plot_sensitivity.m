rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
ROIs = rois;
%% Pull individual ROI parameters and save to CSV and build figure
load(fullfile(datafolder,'avg_hrffits.mat'));
load(fullfile(datafolder,'avg_hrf_roidelta.mat'));
%%
% process con/coh responses
x = 0:.01:1;
clear conRmax lincoh
for si = 1:length(nSIDs)
    for ri = 1:8
        conRmax(si,ri) = sfits{si}{1,4}.roifit{ri}.params.conRmax;
        conc50(si,ri) = sfits{si}{1,4}.roifit{ri}.params.conc50; 
        lincoh(si,ri) = x'\cohModel(x,sfits{si}{1,4}.roifit{ri}.params)';
        onset(si,ri) = sfits{si}{1,4}.roifit{ri}.params.offset;
    end
end

clear data hrfs
params = {'tau1','amp2','tau2'};
for si = 1:length(nSIDs)
    for ri = 1:8
        for pi = 1:length(params)
            data(si,ri,pi) = afits{si}{ri}.params.(params{pi});
        end
        hrfs(si,ri,:) = cc_gamma(0.25:.5:40.5,afits{si}{ri}.params);
    end
end

figure;
for param = 1:3
    subplot(3,1,param)
    dat = data(:,:,param);
    hist(dat(:));
end

%% compute the delta parameter
% for si = 1:length(nSIDs)
%     for ri = 1:8
%         delta(si,ri) = afits{si}{ri}.params.hrfexp;
%     end
% end

d = 2.^cdeltas;

h = figure; hold on

% Generate bar graph style for the parameter estimates
for ri = 1:8
    plot(repmat(ri,1,size(d,1)),d(:,ri),'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerSize',2);
%     mcon = squeeze(mean(bootci(1000,@mean,d(:,ri))));
    plot(ri,mean(d(:,ri)),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);
end
set(gca,'XTick',1:8,'XTickLabel',ROIs,'YTick',[ 0.5 1],'YTickLabel',[ 50 100]);
ylabel('Proportion of linear response (%)');
axis([1 8 0.5 1.1]);
drawPublishAxis('figSize=[8.9,2.25]');

% savepdf(h,fullfile(datafolder,'avg_fmri','delta_sensitivity.pdf'));

%% Compute statistics for the deltas
clear ci
for ri = 1:8
    ci(:,ri) = bootci(10000,@mean,d(:,ri));
    disp(sprintf('%s %1.0f%%, 95%% CI [%1.2f, %1.2f]',rois{ri},mean(ci(:,ri))*100,ci(1,ri)*100,ci(2,ri)*100));
end
disp(sprintf('Mean percentage of linear %1.2f%%',mean(d(:))*100));

%% Plots
cmap = brewermap(7,'PuOr');
h = figure; hold on

% Generate bar graph style for the parameter estimates
for ri = 1:8
    plot(repmat(ri,1,size(conRmax,1)),conRmax(:,ri),'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerSize',2);
    mcon = squeeze(mean(bootci(1000,@mean,conRmax(:,ri))));
    plot(ri,mcon,'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:),'MarkerSize',5);
end
set(gca,'XTick',1:8,'XTickLabel',rois,'YTick',[0 0.25 0.5 1 2]);
ylabel('Sensitivity');
axis([1 8 0 5]);
drawPublishAxis('figSize=[8.9,2.25]');

% savepdf(h,fullfile(datafolder,'avg_fmri','conRmax_sensitivity.pdf'));


h = figure; hold on

% Generate bar graph style for the parameter estimates
for ri = 1:8
    plot(repmat(ri,1,size(conc50,1)),conc50(:,ri),'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerSize',2);
    mcon = squeeze(mean(bootci(1000,@mean,conc50(:,ri))));
    plot(ri,mcon,'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:),'MarkerSize',5);
end
set(gca,'XTick',1:8,'XTickLabel',rois,'YTick',[0 0.25 0.5 1]);
ylabel('Sensitivity');
axis([1 8 0 1]);
drawPublishAxis('figSize=[8.9,2]');

% savepdf(h,fullfile(datafolder,'avg_fmri','conc50_sensitivity.pdf'));


h = figure; hold on

% Generate bar graph style for the parameter estimates
for ri = 1:8
    plot(repmat(ri,1,size(lincoh,1)),lincoh(:,ri),'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerSize',2);
    mcoh = squeeze(mean(bootci(1000,@mean,lincoh(:,ri))));
    plot(ri,mcoh,'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor',cmap(6,:),'MarkerSize',5);
end
set(gca,'XTick',1:8,'XTickLabel',rois,'YTick',[0 0.25 0.5]);
ylabel('Sensitivity');
axis([1 8 0 0.75]);
drawPublishAxis('figSize=[8.9,2.25]');

% savepdf(h,fullfile(datafolder,'avg_fmri','coh_sensitivity.pdf'));

h = figure; hold on

% Generate bar graph style for the parameter estimates
for ri = 1:8
    plot(repmat(ri,1,size(onset,1)),onset(:,ri),'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerSize',2);
    mcoh = squeeze(mean(bootci(1000,@mean,onset(:,ri))));
    plot(ri,mcoh,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',5);
end
set(gca,'XTick',1:8,'XTickLabel',rois,'YTick',[0 0.25 0.5]);
ylabel('Onset parameter');
axis([1 8 0 0.75]);
drawPublishAxis('figSize=[8.9,2]');

% savepdf(h,fullfile(datafolder,'avg_fmri','onset_sensitivity.pdf'));

%% Relative sensitivity

h = figure; hold on

relative = conRmax ./ lincoh;
relative(lincoh<0.01) = NaN;
% Generate bar graph style for the parameter estimates
for ri = 1:8
    plot(repmat(ri,1,size(relative,1)),relative(:,ri),'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerSize',2);
    rel = squeeze(mean(bootci(1000,@nanmean,relative(:,ri))));
    plot(ri,rel,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',5);
end
set(gca,'XTick',1:8,'XTickLabel',rois,'YTick',[0.1 1 10 100 1000]);
set(gca,'YScale','log');
ylabel('Ratio of contrast and coherence parameters');
axis([1 8 0.1 1000]);
drawPublishAxis('figSize=[16,5]');
savepdf(h,fullfile(datafolder,'avg_fmri','relative_sensitivity.pdf'));

%% Stats abpit stiff
clear ci
for ri = 1:8
    ci(:,ri) = bootci(10000,@mean,conc50(:,ri));
    disp(sprintf('%s %1.2f, 95%% CI [%1.2f, %1.2f]',rois{ri},mean(ci(:,ri)),ci(1,ri),ci(2,ri)));
end
disp(sprintf('c50 average %1.2f%%',mean(conc50(:))*100));

%% Onset stats
clear ci
for ri = 1:8
    ci(:,ri) = bootci(10000,@mean,onset(:,ri));
    disp(sprintf('%s %1.2f, 95%% CI [%1.2f, %1.2f]',rois{ri},mean(ci(:,ri)),ci(1,ri),ci(2,ri)));
end
disp(sprintf('c50 average %1.2f%%',mean(onset(:))*100));

%% Covariance matrix of parameters

mat = [d(:) conRmax(:) conc50(:) onset(:) lincoh(:)];

%% 
conRmax_ = conRmax ./ onset;
conRmax_ = conRmax_([1 2 4:11],:);
lincoh_ = lincoh ./ onset;
lincoh_ = lincoh_([1 2 4:11],:);
%% Test of normalization by onset

cmap = brewermap(7,'PuOr');
h = figure; hold on

% Generate bar graph style for the parameter estimates
for ri = 1:8
    plot(repmat(ri,1,size(conRmax_,1)),conRmax_(:,ri),'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerSize',2);
    mcon = squeeze(mean(bootci(1000,@median,conRmax_(:,ri))));
    plot(ri,mcon,'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:),'MarkerSize',5);
end
set(gca,'XTick',1:8,'XTickLabel',rois,'YTick',[0 0.25 0.5 1 2]);
ylabel('Sensitivity');
axis([1 8 0 20]);
drawPublishAxis('figSize=[8.9,2.25]');

% savepdf(h,fullfile(datafolder,'avg_fmri','conRmax_sensitivity.pdf'));

h = figure; hold on

% Generate bar graph style for the parameter estimates
for ri = 1:8
    plot(repmat(ri,1,size(lincoh_,1)),lincoh_(:,ri),'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerSize',2);
    mcoh = squeeze(mean(bootci(1000,@median,lincoh_(:,ri))));
    plot(ri,mcoh,'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor',cmap(6,:),'MarkerSize',5);
end
set(gca,'XTick',1:8,'XTickLabel',rois,'YTick',[0 0.25 0.5]);
ylabel('Sensitivity');
axis([1 8 0 10]);
drawPublishAxis('figSize=[8.9,2.25]');

% savepdf(h,fullfile(datafolder,'avg_fmri','coh_sensitivity.pdf'));
