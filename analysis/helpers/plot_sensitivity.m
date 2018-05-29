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

savepdf(h,fullfile(datafolder,'avg_fmri','conRmax_sensitivity.pdf'));


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
drawPublishAxis('figSize=[8.9,2.25]');

savepdf(h,fullfile(datafolder,'avg_fmri','conc50_sensitivity.pdf'));


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

savepdf(h,fullfile(datafolder,'avg_fmri','coh_sensitivity.pdf'));
