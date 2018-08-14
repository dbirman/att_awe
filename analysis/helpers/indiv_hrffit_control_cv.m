%% Load decon data
% The idea is to do exactly what we were doing with fitCCTimecourseModel
% but just on the deconvolved HRFs themselves and not on the full
% timecourse. This of course has the advantage of fitting infinitely
% better.

% the HRF will be fit as the average of all the 2.5 s conditions (re-scaled to 1)
ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
rOpts = {'resp_prf','resp_all','resp_2','resp_25'};

for si = 1:length(nSIDs)
    subj = nSIDs(si);
    load(fullfile(datafolder,sprintf('s%04.0f_decon.mat',subj)));
    % stored in decondata
    
    % separate the data
    decondata = cvsplit_decondata(decondata);
    
    % for each split
    for fi = 1:decondata.nfolds
        fold = decondata.(sprintf('fold%i',fi));
        % reform correctly
        fold_train = reformAllData_cv(fold,'train');
        fold_test = reformAllData_cv(fold,'test');
        % fit HRF model
        ffit = fitCCHRFModel_hrf(fold_train,'','resp_25');
        % obtain the impulse fit for this training round
        fold_train.hrf = ffit.impulse;
        % full model
        rfit = fitCCHRFModel(fold_train,'fitroi,naka','','resp_25');
        % fit test set
        fold_test.fit = rfit;
        testfit = fitCCHRFModel(fold_test,'predict');
    end
end

%% Prefit the HRFs (Full Fit)
clear hfit hrfs
hrfs = zeros(11,4,81);
dataopts = {'resp_prf','resp_all','resp_2','resp_25'};
parfor si = 1:length(nSIDs)
    cfit = cell(1,4);
    chrfs = zeros(4,81);
    for di = 1:length(dataopts)
        decondata = decondatas{si}.decondata;
        allData = reformAllData(decondata);
        cfit{di} = fitCCHRFModel_hrf(allData,'',dataopts{di});
        chrfs(di,:) = cfit{di}.impulse;
    end
    hrfs(si,:,:) = chrfs;
    hfit{si} = cfit;
end

save(fullfile(datafolder,'avg_hrf.mat'),'hrfs');

%% plot hrfs for comparison
figure;
for si = 1:length(nSIDs)
    subplot(length(nSIDs),1,si); hold on
    load(fullfile(datafolder,'avg_hrf_exp.mat'));
    plot(hrfs(si,:),'-r');
    load(fullfile(datafolder,'avg_hrf.mat'));
    plot(hrfs(si,:),'-b');
end

%%
figure
plot(((.5:1:4).^.01));
%% Fit the exponent?
% % exp = zeros(1,11);
% % parfor si = 1:length(nSIDs)
% %     decondata = decondatas{si}.decondata;
% %     allData = reformAllData(decondata);
% %     allData.hrf = hfit{si}.impulse;
% %     efit{si} = fitCCHRFModel(allData,'fitexp');
% %     exp(si) = hfit{si}.;
% % end
% % 
% % save(fullfile(datafolder,'avg_hrf.mat'),'hrfs');

%% Show hrfs
load(fullfile(datafolder,'avg_hrf.mat'));

h = figure; 

t = 0.25:.5:40.5;
for si = 1:length(nSIDs)
    subplot(length(nSIDs),1,si); hold on
    axis([0 20 -0.5 1]);
%     plot(t,hfit{si}.impulse);
    plot(t,hrfs(si,:),'r');
end

%% RUN
warning('Add cross-validation step to model fitting');
sfits = cell(1,length(nSIDs));
models = {'fitroi,naka','fitroi,linear'};
% models = {'fitroi,linear','fitroi,linear,interaction'};

parfor si = 1:length(nSIDs)
    
    decondata = decondatas{si}.decondata;
    allData = reformAllData(decondata);
    
    % other models: ,'fitroi,naka,nooffset','fitroi,linear,nooffset''fitroi,naka,doubleoffset','fitroi,linear,doubleoffset'
    
    mfit = cell(size(models));
    
    sBICs = zeros(length(models),8);
    for mi = 1:length(models)
        mfit{mi} = fitCCHRFModel(allData,models{mi});
        sBICs(mi,:) = mfit{mi}.AIC;
    end
    
    % find the best NAKA model (use average diff acroiss ROIs)
%     if mean(sBICs(1,:)-sBICs(2,:))<0 % regular better
%         mfit{end+1} = mfit{1};
%     else
%         mfit{end+1} = mfit{2};
%     end
%     % find the best LINEAR model (use average diff acroiss ROIs)
%     if mean(sBICs(3,:)-sBICs(4,:))<0 % regular better
%         mfit{end+1} = mfit{3};
%     else
%         mfit{end+1} = mfit{4};
%     end
%     
    sfits{si} = mfit;
end

% models{end+1} = 'best,naka';
% models{end+1} = 'best,linear';

save(fullfile(datafolder,'avg_hrffits.mat'),'sfits');

%% Interaction test
clear AICs ib
for si = 1:length(nSIDs)
    for mi = 1:2
        for ri = 1:8
            AICs(si,mi,ri) = sfits{si}{mi}.AIC(ri);
            ib(si,ri) = sfits{si}{2}.roifit{ri}.params.inbeta;
        end
    end
end
AICs = squeeze(AICs(:,1,:)-AICs(:,2,:)); % positive = evidence for INTERACTION

%% Load 
load(fullfile(datafolder,'avg_hrffits.mat'));
%% Convert BIC

BICs = zeros(11,length(models),8);
AICs = zeros(11,length(models),8);
for si = 1:length(nSIDs)
    for mi = 1:length(models)
        for ri = 1:8
            BICs(si,mi,ri) = sfits{si}{mi}.BIC(ri);
            AICs(si,mi,ri) = sfits{si}{mi}.AIC(ri);
        end
    end
end

BICs = BICs - repmat(BICs(:,1,:),1,length(models),1);
AICs = AICs - repmat(AICs(:,1,:),1,length(models),1);

AICs_ = squeeze(bootci(1000,@mean,AICs));
mAICs = squeeze(mean(AICs_));
sAICs = squeeze(AICs_(2,:,:))-mAICs;

h = figure; hold on
% subplot(211)
cmap = brewermap(length(models),'Dark2');
clear p
for mi = 1:length(models)
    errbar((1:8)+mi*0.05,mAICs(mi,:),sAICs(mi,:),'-','Color',cmap(mi,:));
    p(mi) = plot((1:8)+mi*0.05,mAICs(mi,:),'ok','MarkerFaceColor',cmap(mi,:),'MarkerEdgeColor','w');
end

% set(gca,'yscale','log');

set(gca,'XTick',1:8,'XTickLabel',ROIs);
% subplot(212)
% plot(BICs');
l = legend(p,models,'Location','NorthEast');
set(l,'box','off','Xcolor',[1 1 1],'Ycolor',[1 1 1]);
 
drawPublishAxis('figSize=[8.5,6]');
savepdf(h,fullfile(datafolder,'avg_models','fmri_model_comp.pdf'));

%% Average R^2 improvement
for si = 1:length(nSIDs)
    r2(si,:) = sfits{si}{2}.r2-sfits{si}{1}.r2;
end

ci = bootci(10000,@mean,r2);
mean(ci)
mean(ci(:))

%% Bar plot for MT
models = {'naka','lin'};
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
BICs = zeros(11,length(models),8);
AICs = zeros(11,length(models),8);
for si = 1:length(nSIDs)
    for mi = 1:length(models)
        for ri = 1:8
            BICs(si,mi,ri) = sfits{si}{mi}.BIC(ri);
            AICs(si,mi,ri) = sfits{si}{mi}.AIC(ri);
        end
    end
end
% use only V3a and MT
clear diffs
AICs = AICs(:,:,[5 8]);
lrois = rois([5 8]);
for ri = 1:2
    lA = AICs(:,:,ri);
    diffs(:,ri) = lA(:,1)-lA(:,2); % if diff > 0 then linear is the better model 
end

% upper plot will be a horizontal bar plot of the # of subjects preferring
% X vs Y at size 4.45 / 2
h = figure;
diffs_v3a = diffs(:,1);
naka = sum(diffs_v3a<-2);
lin = sum(diffs_v3a>2);
none = length(diffs_v3a)-naka-lin;
rectangle('Position',[-2 ri-0.25 4 0.5],'FaceColor',[0.5 0.5 0.5]);
text(0.75,ri,num2str(none),'HorizontalAlignment','center');
if lin>0
    rectangle('Position',[2 ri-0.25 lin 0.5],'FaceColor',[0.75 0.75 0.75]);
    text(2.5,ri,num2str(lin),'HorizontalAlignment','center');
end
if naka>0
    rectangle('Position',[(-2-naka) ri-0.25 naka 0.5],'FaceColor',[0.75 0.75 0.75]);
    text(-2.5,ri,num2str(naka),'HorizontalAlignment','center');
end
%     ps(2) = rectangle('Position',[xc, -2, 0.5, 4],'FaceColor',cmap(2,:));
%     text(xc+.25,0.75,num2str(AIC_o(ci,2)),'HorizontalAlignment','center');

vline(0,'--r');
axis([-9 9 1 3]);
a = axis;
set(gca,'XTick',[-5 5],'XTickLabel',{'Naka-Rushton','Linear'});
set(gca,'YTick',[0],'YTickLabeL',{''});
% set(gca,'XTick',1:2,'XTickLabel',lrois);
% set(gca,'YTick',[-5 5],'YTickLabel',{'Naka-Rushton','Linear'});
xlabel('Evidence for model (\Delta AIC)');
title('Cortical area V3a');
drawPublishAxis('figSize=[4.45,2]');

% savepdf(h,fullfile('~/proj/att_awe/talks/data_figures/v3a_mt_fmri.pdf'));
savepdf(h,fullfile(datafolder,'avg_models','v3a_fmri_bar.pdf'));

% histogram plot
h = figure; hold on

[b,i] = hist(r2(:,5)*100,20);
bins = unique(round(i*10)/10); db = diff(bins);
bins = [bins(1)-db(1)*1.5 bins-db(1)/2];
[b,i] = hist(r2(:,5)*100,bins);
d = diff(i); d = d(1);
for ii = 1:length(b)
    if b(ii)>0
        if abs(i(ii))>2
            rectangle('Position',[i(ii)-d/2 0 d b(ii)],'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0]);
        else
            rectangle('Position',[i(ii)-d/2 0 d b(ii)],'FaceColor',[0.5 .5 .5],'EdgeColor',[0 0 0]);
        end
    end
end
vline(0,'--r');
% a2 = axis;
% mult = max(a2(1)/a(1));
% set(gca,'XTick',[-5 5]*mult,'XTickLabel',{'Naka-Rushton','Linear'});
% set(gca,'XTick',round([-5 5]*mult));
set(gca,'XTick',[-1 0 1]);
xlabel('\Delta R^2 (Linear - Naka-Rushton, %)');
set(gca,'YTick',[0 1 2]);
ylabel('Observers');
axis([a(1)*mult a(2)*mult a2(3) a2(4)]);
drawPublishAxis('figSize=[4.45,2]');
savepdf(h,fullfile(datafolder,'avg_models','v3a_fmri_hist.pdf'));

h = figure;
diffs_mt = diffs(:,2);
naka = sum(diffs_mt<-2);
lin = sum(diffs_mt>2);
none = length(diffs_mt)-naka-lin;
rectangle('Position',[-2 ri-0.25 4 0.5],'FaceColor',[0.5 0.5 0.5]);
text(0.75,ri,num2str(none),'HorizontalAlignment','center');
if lin>0
    rectangle('Position',[2 ri-0.25 lin 0.5],'FaceColor',[0.75 0.75 0.75]);
    text(2.5,ri,num2str(lin),'HorizontalAlignment','center');
end
if naka>0
    rectangle('Position',[(-2-naka) ri-0.25 naka 0.5],'FaceColor',[0.75 0.75 0.75]);
    text(-2.5,ri,num2str(naka),'HorizontalAlignment','center');
end
%     ps(2) = rectangle('Position',[xc, -2, 0.5, 4],'FaceColor',cmap(2,:));
%     text(xc+.25,0.75,num2str(AIC_o(ci,2)),'HorizontalAlignment','center');

vline(0,'--r');
axis([-10 6 1 3]);
a = axis;
set(gca,'XTick',[-5 5],'XTickLabel',{'Naka-Rushton','Linear'});
set(gca,'YTick',[0],'YTickLabeL',{''});
% set(gca,'XTick',1:2,'XTickLabel',lrois);
% set(gca,'YTick',[-5 5],'YTickLabel',{'Naka-Rushton','Linear'});
xlabel('Evidence for model (\Delta AIC)');
title('Cortical area MT');
drawPublishAxis('figSize=[4.45,2]');
savepdf(h,fullfile(datafolder,'avg_models','v3a_mt_bar.pdf'));


% lower plot will be a histogram of the actual AIC values across subjects
% all plots render separately (4.45/2 size) because of mis-aligned axes

% histogram plot
h = figure; hold on

[b,i] = hist(r2(:,8)*100,20);
bins = unique(round(i*10)/10); db = diff(bins);
bins = [bins(1)-db(1)*1.5 bins-db(1)/2];
[b,i] = hist(r2(:,8)*100,bins);
d = diff(i); d = d(1);
for ii = 1:length(b)
    if b(ii)>0
        if abs(i(ii))>2
            rectangle('Position',[i(ii)-d/2 0 d b(ii)],'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0]);
        else
            rectangle('Position',[i(ii)-d/2 0 d b(ii)],'FaceColor',[0.5 .5 .5],'EdgeColor',[0 0 0]);
        end
    end
end
vline(0,'--r');
% a2 = axis;
% mult = max(a2(1)/a(1));
% set(gca,'XTick',[-5 5]*mult);
set(gca,'XTick',[-1 0 1],'XTickLabel',[-1 0 1]);
xlabel('\Delta R^2 (Linear - Naka-Rushton, %)');
set(gca,'YTick',[0 1 2]);
ylabel('Observers');
% axis([a(1)*mult a(2)*mult a2(3) a2(4)]);
drawPublishAxis('figSize=[4.45,2]');
savepdf(h,fullfile(datafolder,'avg_models','v3a_mt_hist.pdf'));

%% Get interaction parameter
inp = zeros(11,length(models),8);
for si = 1:length(nSIDs)
    for mi = 1:length(models)
        for ri = 1:8
            inp(si,mi,ri) = sfits{si}{mi}.roifit{ri}.params.inbeta;
        end
    end
end