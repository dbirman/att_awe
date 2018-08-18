nSIDs = [305 329 43 25 300 346 343 338 344 348];
ncorrespond = [1 2 3 4 5 6 7 9 10 11];

%% Load

fname = fullfile(datafolder,'avg_deconatt','avg_deconEffects_att_cross.mat');
load(fname);

% avgdecon is now available for use

%% Convert allCon to betas

allResp = avgdecon.allResp;
beta = zeros(size(allResp,1),size(allResp,2),size(allResp,3));

load(fullfile(datafolder,'avg_hrf.mat'));

for ni = 1:10
    for ri = 1:8
        for ci = 1:32
            % get the right hrf, i.e. the average of the other observers
            pull = setdiff(1:11,ncorrespond(ni));
            ahrf = squeeze(mean(hrfs(pull,4,:)));
            % ahrf = mean(hrfs,1);
            ahrf = ahrf./max(ahrf);

            beta(ni,ri,ci) = ahrf\squeeze(allResp(ni,ri,ci,:));
        end
    end
end

%% Plot Betas (to look at response functions across conditions)
% code is below under CONTRAST/COHERENCE
conidx = avgdecon.conidx;
cohidx = avgdecon.cohidx;
contrasts = avgdecon.contrasts;
coherences = avgdecon.coherences;
taskidx = avgdecon.taskidx;

%% Fit baseline shift / gain models
% In betaCon and betaCoh we have the fits of the standard HRF (ahrf, loaded
% from avg_hrf.mat). We now want to take these values and basically compare
% them to the values that were saved from the fixation data, and see what
% offset and gain parameters we need to adjust. We'll use the function
% fitCCHRFModel_att (which is dervied from fitCCHRFModel a combination of
% fitCCTimecourseModel). We can do a BIC comparison as well to make sure
% we're fitting the best possible model but it doesn't really matter in
% practice. The end result of this is that there are no differences in the
% attention conditions (excepted possible baseline shifts). We can then
% draw up response function plots comparing, for each ROI,
% attend/unattended/fixated response functions--these can then be compared
% against the same functions generated the behavioral data under various
% assumptions (baseline shift, response gain, or noise reduction or
% something). We find that none of the models fit well suggesting that
% feature-based attention doesn't act as sensory enhancement (in either a 
% bias or gain dependent manner) 

avgdecon.beta = beta;
avgdecon.ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

attfits = cell(1,length(nSIDs));
% fit only the full model
models = {'doublebaseline'};
parfor ni = 1:length(nSIDs)
    % compute fit values for each subject
    ldecon = avgdecon;
    ldecon.beta = beta(ni,:,:);
    
    % try three different fits
    fits = cell(size(models));
    for mi = 1:length(models)
        fits{mi} = fitCCHRFModel_att_cross(ldecon,nSIDs(ni),models{mi},1); % fit contrast gains
    end
    
    
    attfits{ni} = fits;
end

save(fullfile(datafolder,'avg_att_cross_fits.mat'),'attfits');

%% One baseline model
models = {''};
parfor ni = 1:length(nSIDs)
    % compute fit values for each subject
    ldecon = avgdecon;
    ldecon.beta = beta(ni,:,:);
    
    % try three different fits
    fits = cell(size(models));
    for mi = 1:length(models)
        fits{mi} = fitCCHRFModel_att_cross(ldecon,nSIDs(ni),models{mi},1); % fit contrast gains
    end
    
    
    attfits{ni} = fits;
end

save(fullfile(datafolder,'avg_att_cross_fits_sb.mat'),'attfits');

% for ni = 1:length(nSIDs)
%     fits = attfits{ni};
%     save(fullfile(datafolder,sprintf('s%04.0f_fithrf_att.mat',nSIDs(ni))),'fits');
% end

%%
load(fullfile(datafolder,'avg_att_cross_fits.mat'));

%% Bootstrap slopes

% baseline
ocon = zeros(length(nSIDs),8);
ocoh = zeros(length(nSIDs),8);
for ni = 1:length(nSIDs)
    for ri = 1:8
        ocon(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.offset_con;
        ocoh(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.offset_coh;
    end
end

cc_slope('Baseline',ocon,ocoh);

% con
ocon = zeros(length(nSIDs),8);
ocoh = zeros(length(nSIDs),8);
for ni = 1:length(nSIDs)
    for ri = 1:8
        ocon(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.con_congain;
        ocoh(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.con_cohgain;
    end
end

cc_slope('Contrast gain',ocon,ocoh);

% coh
ocon = zeros(length(nSIDs),8);
ocoh = zeros(length(nSIDs),8);
for ni = 1:length(nSIDs)
    for ri = 1:8
        ocon(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.coh_congain;
        ocoh(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.coh_cohgain;
    end
end

cc_slope('Coherence gain',ocon,ocoh);


% I have no idea what the next  400 lines of code do.

%% Check for correlations (baseline)
% ocon = zeros(length(nSIDs),8);
% ocoh = zeros(length(nSIDs),8);
% for ni = 1:length(nSIDs)
%     for ri = 1:8
%         ocon(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.offset_con;
%         ocoh(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.offset_coh;
%     end
% end
% 
% cor = corr(ocon(:),ocoh(:));
% 
% % fit model:
% y = ocoh(:);
% x = ocon(:); x = [ones(size(x)) x];
% 
% b = x\y;
% figure; hold on
% plot(ocon(:),ocoh(:),'*');
% xs = [-0.5 1.5];
% plot(xs,b(1)+xs*b(2),'-r');
% 
% %% Contrast gain
% 
% ocon = zeros(length(nSIDs),8);
% ocoh = zeros(length(nSIDs),8);
% for ni = 1:length(nSIDs)
%     for ri = 1:8
%         ocon(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.con_congain;
%         ocoh(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.con_cohgain;
%     end
% end
% % fit model:
% y = ocoh(:);
% x = ocon(:); x = [ones(size(x)) x];
% 
% b = x\y;
% figure; hold on
% plot(ocon(:),ocoh(:),'*');
% axis([-5 5 -5 5]);
% title('Contrast Gain');
% xlabel('Contrast discrimination')
% ylabel('Coherence discrimination');
% xs = [-0.5 1.5];
% plot(xs,b(1)+xs*b(2),'-r');
% plot([-5 5],[-5 5],'--k');
% 
% %% Coherence gain
% 
% ocon = zeros(length(nSIDs),8);
% ocoh = zeros(length(nSIDs),8);
% for ni = 1:length(nSIDs)
%     for ri = 1:8
%         ocon(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.coh_congain;
%         ocoh(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.coh_cohgain;
%     end
% end
% % fit model:
% y = ocoh(:);
% x = ocon(:); x = [ones(size(x)) x];
% 
% 
% b = x\y;
% figure; hold on
% plot(ocon(:),ocoh(:),'*');
% 
% axis([-5 5 -5 5]);
% xs = [-0.5 1.5];
% title('Coherence Gain');
% xlabel('Contrast discrimination')
% ylabel('Coherence discrimination');
% plot(xs,b(1)+xs*b(2),'-r');
% plot([-5 5],[-5 5],'--k');
% 
% %% Attention model figures: compare offset in acon vs acoh
% 
% load(fullfile(datafolder,'avg_att_fits.mat'));
% ocon = zeros(length(nSIDs),8);
% ocoh = zeros(length(nSIDs),8);
% for ni = 1:length(nSIDs)
%     for ri = 1:8
%         ocon(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.offset_con;
%         ocoh(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.offset_coh;
%     end
% end
% 
% %% Baseline stats
% 
% ci = bootci(1000,@mean,ocon(:,[1 8]));
% mi = bootci(1000,@mean,ocoh(:,[1 8]));
% 
% for ni = 1:10
%     [cr(ni), cp(ni)] = corr(ocon(ni,:)',ocoh(ni,:)');
% end
% 
% %% Baseline plot
% 
% % Test: is there are difference in baseline offset between the attend contrast and
% % attend coherence conditions
% 
% % get info
% bins = floor(4*min([ocon(:) ; ocoh(:)]))/4:.25:ceil(4*max([ocon(:) ; ocoh(:)]))/4;
% 
% h = figure; hold on
% 
% % subplot(121); hold on
% subplot(3,3,[4 5 7 8]); hold on
% plot(ocon(:),ocoh(:),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
% % ,'o',);
% y = ocoh(:);
% x = ocon(:); x = [ones(size(x)) x];
% b = x\y;
% xs = [-0.5 1.5];
% plot([-1 2],[-1 2],'--k');
% plot(xs,b(1)+b(2)*xs,'-r');
% 
% % p = polyfit(ocon(:),ocoh(:),1);
% % x = [min(ocon(:)) max(ocon(:))];
% % y = p(1)*x + p(2);
% % plot(x,y,'--r');
% axis([-1 2 -1 2]);
% % axis equal
% hline(0,'--k');
% vline(0,'--k');
% cc = corrcoef(ocon(:),ocoh(:));
% xlabel('Discriminate contrast');
% ylabel('Discriminate motion');
% 
% drawPublishAxis;
% 
% subplot(3,3,[1 2]); hold on
% title(sprintf('Baseline shift r = %0.2f',cc(1,2)));
% [c,b] = hist(ocon(:),bins);
% bar(b,c,'k');
% vline(0,'--k');
% a = axis;
% axis([-1 2 a(3) a(4)]);
% % set(gca,'XTick',[],'YTick',[]);
% 
% drawPublishAxis;
% 
% subplot(3,3,[6 9]); hold on
% [c,b] = hist(ocoh(:),bins);
% barh(b,c,'k');
% hline(0,'--k');
% a = axis;
% axis([a(1) a(2) -1 2]);
% % drawPublishAxis;
% 
% drawPublishAxis('figSize=[4.5 4.5]');
% 
% savepdf(h,fullfile(datafolder,'all_att','baseline.pdf'));
% %% Compare congain in acon vs. acoh
% ocon = zeros(length(nSIDs),8);
% ocoh = zeros(length(nSIDs),8);
% for ni = 1:length(nSIDs)
%     for ri = 1:8
%         ocon(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.con_congain;
%         ocoh(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.con_cohgain;
%     end
% end
% bins = floor(4*min([ocon(:) ; ocoh(:)]))/4:.8:ceil(4*max([ocon(:) ; ocoh(:)]))/4;
% 
% 
% % Test: is there are difference in baseline offset between the attend contrast and
% % attend coherence conditions
% 
% h = figure;
% subplot(3,3,[4 5 7 8]); hold on
% plot(ocon(:),ocoh(:),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
% hline(1,'--k');
% vline(1,'--k');
% % p = polyfit(ocon(:),ocoh(:),1);
% % x = [min(ocon(:)) max(ocon(:))];
% % y = p(1)*x + p(2);
% % plot(x,y,'--r');
% 
% y = ocoh(:);
% x = ocon(:); x = [ones(size(x)) x];
% b = x\y;
% xs = [-3 3];
% plot([-5 5],[-5 5],'--k');
% plot(xs,b(1)+b(2)*xs,'--r');
% 
% cc = corrcoef(ocon(:),ocoh(:));
% xlabel('Discriminate contrast');
% ylabel('Discriminate motion');
% axis([-5 5 -5 5]);
% drawPublishAxis;
% 
% subplot(3,3,[1 2]); hold on
% title(sprintf('Contrast gain r = %0.2f',cc(1,2)));
% [c,b] = hist(ocon(:),bins);
% bar(b,c,'k');
% vline(1,'--k');
% a = axis;
% axis([-5 5 a(3) a(4)]);
% % set(gca,'XTick',[],'YTick',[]);
% 
% drawPublishAxis;
% 
% subplot(3,3,[6 9]); hold on
% [c,b] = hist(ocoh(:),bins);
% barh(b,c,'k');
% hline(1,'--k');
% a = axis;
% axis([a(1) a(2) -5 5]);
% drawPublishAxis('figSize=[4.5 4.5]');
% 
% savepdf(h,fullfile(datafolder,'all_att','congain.pdf'));
% 
% %% Compare cohgain
% ocon = zeros(length(nSIDs),8);
% ocoh = zeros(length(nSIDs),8);
% for ni = 1:length(nSIDs)
%     for ri = 1:8
%         ocon(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.coh_congain;
%         ocoh(ni,ri) = attfits{ni}{1}.roifit{ri}.roiparams.coh_cohgain;
%     end
% end
% bins = floor(4*min([ocon(:) ; ocoh(:)]))/4:.8:ceil(4*max([ocon(:) ; ocoh(:)]))/4;
% 
% 
% % Test: is there are difference in baseline offset between the attend contrast and
% % attend coherence conditions
% 
% h = figure;
% subplot(3,3,[4 5 7 8]); hold on
% plot(ocon(:),ocoh(:),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
% hline(1,'--k');
% vline(1,'--k');
% % p = polyfit(ocon(:),ocoh(:),1);
% % x = [min(ocon(:)) max(ocon(:))];
% % y = p(1)*x + p(2);
% y = ocoh(:);
% x = ocon(:); x = [ones(size(x)) x];
% b = x\y;
% xs = [-3 3];
% plot([-5 5],[-5 5],'--k');
% plot(xs,b(1)+b(2)*xs,'--r');
% 
% 
% cc = corrcoef(ocon(:),ocoh(:));
% xlabel('Discriminate contrast');
% ylabel('Discriminate motion');
% axis([-5 5 -5 5]);
% drawPublishAxis;
% % 
% subplot(3,3,[1 2]); hold on
% title(sprintf('Motion gain r = %0.2f',cc(1,2)));
% [c,b] = hist(ocon(:),bins);
% bar(b,c,'k');
% vline(1,'--k');
% a = axis;
% axis([-5 5 a(3) a(4)]);
% % set(gca,'XTick',[],'YTick',[]);
% 
% drawPublishAxis;
% % 
% subplot(3,3,[6 9]); hold on
% [c,b] = hist(ocoh(:),bins);
% barh(b,c,'k');
% hline(1,'--k');
% a = axis;
% axis([a(1) a(2) -5 5]);
% drawPublishAxis('figSize=[4.5 4.5]');
% % 
% savepdf(h,fullfile(datafolder,'all_att','cohgain.pdf'));
% 
% %% Get all BIC values
% BICs = BICs - repmat(BICs(:,1,:),1,5,1);
% 
% [BICs_min, BICs_mini] = min(BICs,[],2);
% 
% % find best model across subjects (per ROI)
% bestmodels = squeeze(mode(BICs_mini,1));
% 
% %% Plot baseline / gain model effects
% 
% % First convert into a useful space
% % 1 = fixation
% % 2 = resp about motion
% % 3 = resp about contrast
% resp = zeros(length(nSIDs),2,8,3,1001); 
% 
% for ni = 1:length(nSIDs)
%     for ri = 1:8
%         fit = attfits{ni}{1};
%         resp(ni,2,ri,1,:) = fit.roifit{ri}.conresp;
%         resp(ni,2,ri,2,:) = fit.roifit{ri}.conresp_coh;
%         resp(ni,2,ri,3,:) = fit.roifit{ri}.conresp_con;
%         resp(ni,1,ri,1,:) = fit.roifit{ri}.cohresp;
%         resp(ni,1,ri,2,:) = fit.roifit{ri}.cohresp_coh;
%         resp(ni,1,ri,3,:) = fit.roifit{ri}.cohresp_con;
%     end
% end
% 
% % Average across subjects
% resp = squeeze(bootci(1000,@mean,resp));
% mresp = squeeze(mean(resp));
% 
% %% Get ready to plot real data
% 
% conidx = avgdecon.conidx;
% cohidx = avgdecon.cohidx;
% contrasts = avgdecon.contrasts;
% coherences = avgdecon.coherences;
% taskidx = avgdecon.taskidx;
% deltaidx = avgdecon.deltaidx;
% 
% con = conidx(deltaidx==0);
% coh = cohidx(deltaidx==0);
% task = taskidx(deltaidx==0);
% bCon = betaCon(:,:,deltaidx==0);
% bCoh = betaCoh(:,:,deltaidx==0);
% 
% % average
% bCon_ = squeeze(mean(bCon));
% bCoh_ = squeeze(mean(bCoh));
% %% Plot
% x = 0:.001:1;
% h = figure; hold on
% rois = [1 2 3 5 8];
% for rii = 1:length(rois)
%     ri = rois(rii);
%     cmap = brewermap(3,'PuOr');
%     
%     subplot(2,5,rii); hold on
%     % plot fixation as a solid line
%     plot(x,squeeze(mresp(2,ri,1,:)),'-','Color',cmap(1,:),'LineWidth',2);
%     % plot attend contrast as a thin solid
%     plot(x,squeeze(mresp(2,ri,3,:)),'-','Color',cmap(1,:),'LineWidth',1);
%     % plot attend coherence as a thin dashed
%     plot(x,squeeze(mresp(2,ri,2,:)),'--','Color',cmap(1,:),'LineWidth',1);
%     
%     % plot attend data (circles for contrast, triangles for coherence)
%     plot(con(task==2),bCon_(ri,task==2),'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','w');
%     plot(con(task==1),bCon_(ri,task==1),'^','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','w');
%     
%     title(sprintf('%s',ROIs{ri}));
%     
%     if ri==1
%         axis([0 1 0 2.5]);
%         xlabel('\Delta stimulus (%)');
%         ylabel('\Delta signal (%)');
%         set(gca,'XTick',[0 1],'XTickLabel',{'0%','100%'},'YTick',[0 1 2]);
%     else
%         axis([0 1 0 1.5]);
%         set(gca,'XTick',[0 1],'XTickLabel',{'0%','100%'},'YTick',[0 1]);
%     end
%     
%    
%     clear p
%     if ri==5
%         p(1) = plot(-1,-1,'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','w');
%         p(2) = plot(-1,-1,'o','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','w');
%         legend(p,{'Contrast','Coherence'},'FontSize',7,'FontName','Helvetica');
%         l = legend(gca,'boxoff');
%         set(l,'Color','none');
%     end
%     
%     drawPublishAxis;
%     
%     subplot(2,5,rii+5); hold on
%     % plot fixation as a solid line
%     plot(x,squeeze(mresp(1,ri,1,:)),'-','Color',cmap(3,:),'LineWidth',2);
%     % plot attend contrast as a thin solid
%     plot(x,squeeze(mresp(1,ri,3,:)),'-','Color',cmap(3,:),'LineWidth',1);
%     % plot attend coherence as a thin dashed
%     plot(x,squeeze(mresp(1,ri,2,:)),'--','Color',cmap(3,:),'LineWidth',1);
%     
%     % plot attend data (circles for contrast, triangles for coherence)
%     plot(coh(task==2),bCoh_(ri,task==2),'o','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','w');
%     plot(coh(task==1),bCoh_(ri,task==1),'^','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','w');
%     
%     clear p
%     if ri==5
%         p(1) = plot(-1,-1,'-k','LineWidth',5); p(2) = plot(-1,-1,'-ok','LineWidth',1); p(3) = plot(-1,-1,'--^k','LineWidth',1);
% 
%         legend(p,{'Fixate','Respond to Contrast','Respond to Motion'},'FontSize',7,'FontName','Helvetica'); 
%         l = legend(gca,'boxoff');
%         set(l,'Color','none');
%     end
%     
%     if ri==1
%         axis([0 1 0 2.5]);
%         xlabel('\Delta stimulus (%)');
%         ylabel('\Delta signal (%)');
%         set(gca,'XTick',[0 1],'XTickLabel',{'0%','100%'},'YTick',[0 1 2]);
%     else
%         axis([0 1 0 1.5]);
%         set(gca,'XTick',[0 1],'XTickLabel',{'0%','100%'},'YTick',[0 1]);
%     end
%     
%     title(sprintf('%s',ROIs{ri}));
%     
%     
%     drawPublishAxis('figSize=[10.5,5.5]');
% end
% savepdf(h,fullfile(datafolder,'avg_fitatt','avg_resp.pdf'));