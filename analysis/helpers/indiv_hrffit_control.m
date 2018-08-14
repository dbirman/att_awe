%% Run the fitCCHRFModel code on thisc
% The idea is to do exactly what we were doing with fitCCTimecourseModel
% but just on the deconvolved HRFs themselves and not on the full
% timecourse. This of course has the advantage of fitting infinitely
% better.

% the HRF will be fit as the average of all the 2.5 s conditions (re-scaled to 1)
ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

dataopts = {'resp_prf','resp_all','resp_2','resp_25'};

%% Load decondatas for HRF
for si = 1:length(nSIDs)
    subj = nSIDs(si);
    decondatas{si} = load(fullfile(datafolder,sprintf('s%04.0f_decon.mat',subj)));
end

%% Prefit the HRFs (Full Fit)
clear hfit hrfs
hrfs = zeros(11,4,81);
parfor si = 1:length(nSIDs)
    cfit = cell(1,4);
    chrfs = zeros(4,81);
    for di = 1:length(dataopts)
        decondata = decondatas{si}.decondata;
        allData = reformAllData(decondata,si);
        cfit{di} = fitCCHRFModel_hrf(allData,'',dataopts{di});
        chrfs(di,:) = cfit{di}.impulse;
    end
    hrfs(si,:,:) = chrfs;
    hfit{si} = cfit;
end

% save(fullfile(datafolder,'avg_hrf.mat'),'hrfs');
% save(fullfile(datafolder,'avg_hrffits_decon.mat'),'hfit');


%% Prefit the HRFs (Individual ROI fits)
clear hfit hrfs
hrfs = zeros(11,4,81);
parfor si = 1:length(nSIDs)
    cfit = cell(1,4);
    cdelta = zeros(1,8);
    cfits = cell(1,8);
    for ri = 1:8
        decondata = decondatas{si}.decondata;
        allData = reformAllData(decondata,si);
        cfit = fitCCHRFModel_hrf(allData,ROIs{ri},'resp_25');
        cdelta(ri) = cfit.params.hrfexp;
        disp(ri);
        cfits{ri} = cfit;
    end
    cdeltas(si,:) = cdelta;
    afits{si} = cfits;
end

save(fullfile(datafolder,'avg_hrf_roidelta.mat'),'afits','cdeltas');

%% Look at parameters
params = {'hrfexp','tau1','amp2','tau2'};
ps = zeros(4,11);
for si = 1:length(nSIDs)
    for pi = 1:length(params)
        disp(sprintf('%s: %01.3f',params{pi},hfit{si}{4}.params.(params{pi})));
        ps(pi,si) = hfit{si}{4}.params.(params{pi});
    end
end

%% HRF Comparison
hrf_comparison;

%% RUN (restricted to di==4)
sfits = cell(1,length(nSIDs));
models = {'fitroi,exp','fitroi,cohlinear','fitroi,conlinear'};

parfor si = 1:length(nSIDs)
    mfit = cell(length(models),4);

    for di = 1:length(dataopts)        
        decondata = decondatas{si}.decondata;
        allData = reformAllData(decondata,si);

        for mi = 1:length(models)
            mfit{mi,di} = fitCCHRFModel(allData,models{mi},'',dataopts{di},1);
        end
    
    end
    sfits{si} = mfit;
end

% save(fullfile(datafolder,'avg_hrffits_lin4.mat'),'sfits');

save(fullfile(datafolder,'avg_hrffits.mat'),'sfits');

%% Check contrast linearity
load(fullfile(datafolder,'avg_hrffits_lin4.mat'));

for si = 1:length(nSIDs)
    for mi = 1:3
        r2(si,mi,:) = sfits{si}{mi,4}.cv.r2;
    end
end

% compute delta
delta = squeeze(r2(:,3,:)-r2(:,1,:));
ci = bootci(1000,@mean,delta);

ci = ci;
%%
str = '';
for ri = 1:8
    str = strcat(str,sprintf('%s %0.3f, 95%% CI [%0.3f %0.3f]; ',ROIs{ri},mean(ci(:,ri)),ci(1,ri),ci(2,ri)));
end
disp(str);

%% RUN VARIABLE EXPONENT

%% Permutation test
clear pr2
parfor si = 1:length(nSIDs)
    disp(si);
    decondata = decondatas{si}.decondata;
    allData = reformAllData(decondata,si);
    
    for pi = 1:100
        pfit = fitCCHRFModel(allData,'fitroi,exp,permute','','resp_25',1);
    
        pr2(si,pi,:) = pfit.cv.r2; 
        disp(sprintf('%i done %i',si,pi));
    end
end
save(fullfile(datafolder,'permutation_r2.mat'),'pr2');

%% Use permutation results 

load(fullfile(datafolder,'avg_hrffits.mat'));
load(fullfile(datafolder,'permutation_r2.mat'));
%%
for si = 1:11
    for mi = 1:2
        for di = 1:4
            r2(si,mi,di,:) = sfits{si}{mi,di}.cv.r2;
        end
    end
end
        
r2_ = squeeze(mean(r2(:,1,4,:),2));

clear dr2
for pi = 1:100
    dr2(:,:,pi) = r2_-squeeze(pr2(:,pi,:));
end
mean(dr2(:))
dr2 = bootci(1000,@mean,dr2(:));

dr2_ = squeeze(mean(dr2));
dr2_s = squeeze(dr2(2,:,:,:))-dr2_;
%% Load
load(fullfile(datafolder,'avg_hrffits.mat'));

%% Compare cross-val and full r^2 values
clear r2
for ni = 1:length(nSIDs)
    r2(ni,:) = sfits{ni}{1,4}.r2-sfits{ni}{1,4}.cv.r2;
end

%% Plot test data (con/coh: 25/100/250, 75/100/250, 75/100/4000)
% uses time data
clear cresp model
cresp = nan(11,3,81);
model = nan(11,3,81);
for si = 1:11
    time = sfits{si}{1,4}.roifit{1}.cc;
    if size(time.cresp,2)>10
        idxs = [3 12 20];
        cresp(si,:,:) = squeeze(time.cresp(1,idxs,:));
        model(si,:,:) = squeeze(time.model(1,idxs,:));
    end
end
cresp_ = squeeze(nanmean(cresp));
cresp = squeeze(bootci(1000,@nanmean,cresp));
% cresp_ = squeeze(mean(cresp));
cresp_s = squeeze(cresp(2,:,:))-cresp_;
model_ = squeeze(nanmean(model));
model = squeeze(bootci(1000,@nanmean,model));
model_s = squeeze(model(2,:,:))-model_;

h = figure;
for i = 1:3
    subplot(3,1,i); hold on
    plot(0.25:.5:40.5,model_(i,:),'-r');
    plot(0.25:.5:40.5,cresp_(i,:),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
    errbar(0.25:.5:40.5,cresp_(i,:),cresp_s(i,:),'-k');
    axis([0 15 -1 3.5]);
    drawPublishAxis('figSize=[6,10]');
end
savepdf(h,fullfile(datafolder,'model_fig.pdf'));
%% r2 and AIC comparison of models
clear r2
for si = 1:length(nSIDs)
    for mi = 1:length(models)
        for di = 1:length(dataopts)
            r2(si,mi,di,:) = sfits{si}{mi,di}.cv.r2;
            % temp because I fucked up
%             r2(si,mi,di) = myr2(sfits{si}{mi,di}.cv.y,sfits{si}{mi,di}.cv.y_);
        end
    end
end
%% Stats for cohcon paper
ci = bootci(10000,@nanmean,r2);
str = '';
for ri = 1:8
    str = strcat(str,sprintf(' %s, R^2: %0.2f 95%% CI [%0.2f %0.2f];',ROIs{ri},mean(ci(:,ri)),ci(1,ri),ci(2,ri)));
end
disp(str);


%% Linear vs non-linear stats
dexp = squeeze(r2(:,2,4,:)-r2(:,1,4,:));

%% Save r2 into table (for cohcon paper)
csvwriteh(fullfile(datafolder,'hrffit_r2.csv'),round(r2*100)/100,ROIs);

%% plot difference between models
% for mi = 1:4
%     h = figure;
%     r2_ = r2(:,:,mi,5);
%     hold on
%     plot(r2_(:,1),r2_(:,2),'*');
%     plot([min(r2_(:)) max(r2_(:))],[min(r2_(:)) max(r2_(:))],'--k');
%     xlabel('Naka R^2');
%     ylabel('Linear R^2');
%     title(dataopts{mi});
%     naka = sum(r2_(:,1)>r2_(:,2));
%     linear = sum(r2_(:,2)>r2_(:,1));
%     disp(sprintf('In favor of Naka: %i Linear: %i',naka,linear));
%     nci = bootci(10000,@mean,r2_(:,1));
%     nlin = bootci(10000,@mean,r2_(:,2));
%     disp(sprintf('Naka R^2: %01.2f [%01.2f %01.2f]',mean(nci),nci(1),nci(2)));
%     disp(sprintf('Linear R^2: %01.2f [%01.2f %01.2f]',mean(nlin),nlin(1),nlin(2)));
%     dif = r2_(:,1)-r2_(:,2);
%     dci = bootci(10000,@mean,dif);
%     disp(sprintf('Difference Naka-Lin: %01.2f [%01.2f %01.2f]',mean(dci),dci(1),dci(2)));
% end

% %% Figure (using resp-25)
% riopts = [5 8];
% for ri = 1:length(riopts)
%     h = figure;
%     hold on
%     dat = r2(:,:,4,riopts(ri));
%     plot(dat(:,1),dat(:,2),'ok','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
% %     a = axis;
%     plot([0.2 0.6],[0.2 0.6],'--r');
% %     axis(a);
%     axis([0.175 0.625 0.175 0.625]);
%     set(gca,'XTick',[0.2 0.6],'YTick',[0.2 0.6]);
%     xlabel('Naka-Rushton model (R^2)');
%     if ri==1
%         ylabel('Linear model (R^2)');
%     end
%     drawPublishAxis('figSize=[4.45,2]');
%     savepdf(h,fullfile(datafolder,'avg_models',sprintf('%s_cv.pdf',ROIs{riopts(ri)})));
% end

%% Compare the contrast/coherence relative sensitivities across the four voxel selection conditons
voxel_selection_comparison;



%% Compare voxel-choice on con/coh models
cohcon_ratio_comparison;

%% CONTROL EXPERIMENTS

%% Fit interaction model with multiplicative gain

sfits_i = cell(1,length(nSIDs));
parfor si = 1:length(nSIDs)
    % use the voxel-25 model
        
    decondata = decondatas{si}.decondata;
    allData = reformAllData(decondata,si);
    
    sfits_i{si} = fitCCHRFModel(allData,'fitroi,exp,interaction','','resp_25',1);
end
save(fullfile(datafolder,'ifits.mat'),'sfits_i');

%% Fit interaction model with zero condition and non-zero condition, but not just a multiplicative gain
sfits_iz = cell(1,length(nSIDs));
parfor si = 1:length(nSIDs)
    % use the voxel-25 model
        
    decondata = decondatas{si}.decondata;
    allData = reformAllData(decondata,si);
    
    sfits_iz{si} = fitCCHRFModel(allData,'fitroi,exp,intzero','','resp_25',1);
end

save(fullfile(datafolder,'ifits_zero.mat'),'sfits_iz');


%% Interaction test
load(fullfile(datafolder,'ifits_zero.mat'));
for si = 1:length(nSIDs)
%     r2_i(si,:) = sfits_i{si}.cv.r2;
    r2_iz(si,:) = sfits_iz{si}.cv.r2;
end
r2_o = squeeze(r2(:,2,4,:));

clear ci_i
for ri = 1:8
    ci_i(ri,:) = bootci(500,@mean,r2_iz(:,ri)-r2_o(:,ri));
end

out = '';
for ri = 1:8
    out = strcat(out,sprintf(' %s: %01.2f%%,',ROIs{ri},mean(ci_i(ri,:)*100,2)));
end
disp(out);

direction = {'reduces','increases'};
dir = direction{(mean(ci_i)>0)+1};
disp(sprintf('Interaction effect %s R^2 by, on average across ROIs, %01.2f%% [%01.2f %01.2f]',dir,mean(mean(ci_i*100)),mean(ci_i(:,1))*100,mean(ci_i(:,2))*100));

%% Interaction test ignoring subject 3
usids = 1:length(nSIDs);
usids = setdiff(usids,3);

clear r2_i r2_o

for ii = 1:length(usids)
    si = usids(ii);
    r2_i(ii,:) = sfits_iz{si}.cv.r2;
end
r2_o = squeeze(r2(:,2,4,:));
r2_o = r2_o(usids,:);

r2_d = r2_i-r2_o;
r2_d = r2_d*100; % convert to percentage

clear ci_i
for ri = 1:8
    ci_i(ii,:) = bootci(1000,@mean,r2_d(:,ri));
end

out = '';
for ri = 1:8
    out = strcat(out,sprintf(' %s: %01.2f%%,',ROIs{ri},mean(ci_i(ri,:),2)));
end
disp(out);

direction = {'reduces','increases'};
dir = direction{(mean(ci_i)>0)+1};
disp(sprintf('Interaction effect %s R^2 by, on average across ROIs, %01.2f%% [%01.2f %01.2f]',dir,mean(mean(ci_i)),mean(ci_i(:,1)),mean(ci_i(:,2))));

%% Load 
load(fullfile(datafolder,'avg_hrffits.mat'));


%% Get interaction parameter
inp = zeros(11,length(models),8);
for si = 1:length(nSIDs)
    for mi = 1:length(models)
        for ri = 1:8
            inp(si,mi,ri) = sfits_i{si}{mi}.roifit{ri}.params.inbeta;
        end
    end
end

%% Plot the interaction zero functions for comparison
x = 0:.01:1;
clear resp p
for si = 1:length(nSIDs)
    for ri = 1:8
        p{si,ri} = sfits_iz{si}.roifit{ri}.params;
        resp(si,1,ri,1,:) = cohModel(x,p{si,ri});
        resp(si,1,ri,2,:) = conModel(x,p{si,ri});
        r0(si,ri,2) = p{si,ri}.conRmax_0;
        r0(si,ri,1) = p{si,ri}.conRmax;
        % add zero conditions
        p{si,ri}.conc50 = p{si,ri}.conc50_0;
        p{si,ri}.conRmax = p{si,ri}.conRmax_0;
        p{si,ri}.cohalpha = p{si,ri}.cohalpha_0;
        p{si,ri}.cohkappa = p{si,ri}.cohkappa_0;
        resp(si,2,ri,1,:) = cohModel(x,p{si,ri});
        resp(si,2,ri,2,:) = conModel(x,p{si,ri});
        
    end
end

r0 = r0(:,:,2)./r0(:,:,1);
ci = bootci(1000,@mean,r0(:));

resp = resp([1 2 4:(size(resp,1)-1)],:,:,:,:);

resp = squeeze(mean(resp));

h = figure(1);

cmap = brewermap(7,'PuOr');
hold on
clear p
for ri = 1:8
    subplot(4,2,ri); hold on
    p(1) = plot(x,squeeze(resp(1,ri,2,:)),'-','Color',cmap(2,:));
    p(2) = plot(x,squeeze(resp(1,ri,1,:)),'-','Color',cmap(6,:));
    
    p(3) = plot(x,squeeze(resp(2,ri,2,:)),'-','Color',cmap(3,:));
    p(4) = plot(x,squeeze(resp(2,ri,1,:)),'-','Color',cmap(5,:));
    if ri==1
        legend(p,{'Contrast','Coherence','Contrast - no coh','Coherence - no con'});
    end
end


savepdf(h,fullfile(datafolder,'avg_fmri','zerointeraction.pdf'));