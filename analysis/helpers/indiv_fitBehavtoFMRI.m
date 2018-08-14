% Fit the behavioral response models to the fMRI data
%% Load behavioral models

%sids = {'s043','s300','s305','s329','s025'};
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];
sids = {};
for bi = 1:length(aSIDs)
    sids{end+1} = sprintf('s%03.0f',aSIDs(bi));
end
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

% combine data
control = zeros(length(sids),2,4);
x = 0:.001:1;
respfits = zeros(length(sids),2,8,length(x));
for si = 1:length(sids)
    load(fullfile(datafolder,sprintf('%s_data.mat',sids{si})));
    
    % for the optimal behavior we need the response functions from each set
    % of fits data.fits{1:4}. Get them and extrapolate
    
    for i = 1:8
        respfits(si,1,i,:) = cohModel(x,data.fits{i}.params);
        respfits(si,2,i,:) = conModel(x,data.fits{i}.params);
    end
end
%% Average models across subjects

respfits_ = squeeze(bootci(1000,@nanmean,respfits));
respfitz = squeeze(mean(respfits_));

%% Fit these behavioral models to the HRF data
% get subjects and subject data
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];

for si = 1:length(nSIDs)
    subj = nSIDs(si);
    decondatas{si} = load(fullfile(datafolder,sprintf('s%04.0f_decon.mat',subj)));
end
params = [4,2,6,4];
    model_subset = 1:2:7;

sfits = cell(1,length(nSIDs));
for si = 1:length(nSIDs)
    decondata = decondatas{si}.decondata;
    allData = reformAllData(decondata);
        
    
    mfit = cell(1,length(model_subset));
    
    for mi = 1:length(model_subset)
        pfit = struct;
        pfit.x = x;
        pfit.con = squeeze(respfitz(2,model_subset(mi),:));
        pfit.coh = squeeze(respfitz(1,model_subset(mi),:));
        
        mfit{mi} = fitCCHRFModel(allData,'fitsigma,doubleoffset',pfit);
    end
    
    sfits{si} = mfit;
    
end

% save(fullfile(datafolder,'avg_hrffits.mat'),'sfits');

%% Copy over sigmas
sigmas = zeros(11,4,8,2);
BICs = zeros(11,4,8);
likes = zeros(11,4,8);
cc = zeros(11,4,8,20);
time = zeros(11,4,8,20);
x = 0:.001:1;
mresp = zeros(11,4,8,2,length(x));
offsets = zeros(11,4,8,2);
for ni = 1:length(nSIDs)
    for mi = 1:4
        for ri = 1:8
            cp = sfits{ni}{mi}.roifit{ri}.params;
            sigmas(ni,mi,ri,1) = cp.sigmacoh;
            sigmas(ni,mi,ri,2) = cp.sigmacon;
            BICs(ni,mi,ri) = sfits{ni}{mi}.roifit{ri}.BIC;
            likes(ni,mi,ri) = sfits{ni}{mi}.roifit{ri}.like;
            cc(ni,mi,ri,:) = sfits{ni}{mi}.roifit{ri}.cc.resp_;
            if length(sfits{ni}{mi}.roifit{ri}.time.resp_)==20
                time(ni,mi,ri,:) = sfits{ni}{mi}.roifit{ri}.time.resp_;
            end
            mresp(ni,mi,ri,2,:) = squeeze(respfitz(2,model_subset(mi),:));%cohModel(x,cp);
            mresp(ni,mi,ri,1,:) = squeeze(respfitz(1,model_subset(mi),:));
            if isfield(cp,'offset')
                offsets(ni,mi,ri,1) = cp.offset;
                offsets(ni,mi,ri,2) = cp.offset;
            else
                offsets(ni,mi,ri,1) = cp.cohoffset;
                offsets(ni,mi,ri,2) = cp.conoffset;
            end
        end
    end
end

%% Average sigma
sigmas_s = squeeze(bootci(1000,@mean,sigmas));
sigmas_ = squeeze(mean(sigmas_s));
%% Try removing the contrast V1 sigma and look for delta
% remove the sigma for contrast v1, check delta against that
sigmas_con_v1 = sigmas - repmat(sigmas(:,:,1,2),1,1,8,2);
s_c_v1s = squeeze(bootci(1000,@mean,sigmas_con_v1));
s_c_v1 = squeeze(mean(s_c_v1s));

sigmas_coh_mt = sigmas - repmat(sigmas(:,:,8,1),1,1,8,2);
s_c_mts = squeeze(bootci(1000,@mean,sigmas_coh_mt));
s_c_mt = squeeze(mean(s_c_mts));

%% Generate plot comparing con/v1 sigma to coherence rois and coh/mt to contrast rois
h = figure(34); clf

subplot(211); hold on

cmap = brewermap(7,'PuOr');
for mi = 4
    errbar(1:8,abs(squeeze(s_c_v1(mi,:,1))),squeeze(squeeze(s_c_v1s(2,mi,:,1))-s_c_v1(mi,:,1)'),'Color',cmap(6,:));
    plot(1:8,abs(squeeze(s_c_v1(mi,:,1))),'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w');
end
axis([1 8 0 0.12]);
set(gca,'XTick',1:8,'XTickLabel',rois,'YTick',[0 0.1],'YTickLabeL',{'0%','10%'});
hline(0,'--r');
title('Contrast sigma from V1');
ylabel('\Delta sigma (%)');
drawPublishAxis('figSize=[8.5,7]');

subplot(212); hold on
for mi = 4
    errbar(1:8,abs(squeeze(s_c_mt(mi,:,2))),squeeze(squeeze(s_c_mts(2,mi,:,2))-s_c_mt(mi,:,2)'),'Color',cmap(2,:));
    plot(1:8,abs(squeeze(s_c_mt(mi,:,2))),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w');
end
axis([1 8 0 0.12]);
title('Coherence sigma from MT');
set(gca,'XTick',1:8,'XTickLabel',rois,'YTick',[0 0.1],'YTickLabeL',{'0%','10%'});
hline(0,'--r');
ylabel('\Delta sigma (%)');
drawPublishAxis('figSize=[8.5,7]');

% savepdf(h,fullfile(datafolder,'avg_models','behavToFMRI_sigmas.pdf'));

%% Check BIC values 

models = {'Exponential','Linear','Naka','Exp-linear'};
% For each ROI, for each feature, compare all the model values to EXP
BICs = squeeze(BICs(:,:,:));

BICs = BICs - repmat(BICs(:,1,:),1,4,1);

b_con = squeeze(mean(bootci(1000,@mean,BICs)));
% b_coh = squeeze(mean(bootci(1000,@mean,BICs_coh)));

h= figure;
subplot(211);
plot(b_con');
title('Contrast');
legend(models);
subplot(212);
plot(b_coh');
title('Coherence');








%%%%%%%%%%%%%%%%%%%
%% TODO: Generate plots showing behavioral fit to fMRI data %%
%%%%%

%% HRF Fit plots

ccon = sfits{1}{1,1}.roifit{1}.cc.con;
ccoh = sfits{1}{1,1}.roifit{1}.cc.coh;

cc_ = squeeze(mean(bootci(1000,@mean,cc)));
mresp_ = squeeze(mean(bootci(1000,@mean,mresp)));

offsets_ = squeeze(mean(bootci(1000,@mean,offsets)));

%% Create plots
models = {'naka','lin','exp','explin'};
ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

conidx = 6:5:20;
contrasts = [.25 .5 .75];
cohidx = 2:5;
coherences = [.25 .5 .75 1];
timidx = 16:20;
times = [.250 .500 1.000 2.000 4.000];

mopts = 1:4;
% mopts = 1;
for mi=1:length(mopts)
    
    h = figure; clf
    cmap = brewermap(3,'PuOr');
    for ri = 1:8
        crespcon_ = squeeze(mresp_(mi,:,2,:))*sigmas_(mi,ri,2)+offsets_(mi,ri,2);
        crespcoh_ = squeeze(mresp_(mi,:,1,:))*sigmas_(mi,ri,1)+offsets_(mi,ri,1);
        % CONTRAST
        subplot(4,2,ri); hold on
        plot(x,crespcon_(ri,:),'-','Color',cmap(1,:));
%         errbar(contrasts,squeeze(cc_(mi,ri,2,conidx)),squeeze(cc(2,ri,conidx))-cc_(ri,conidx)','-','Color',cmap(1,:));
        plot(contrasts,squeeze(cc_(mi,ri,conidx)),'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','white','MarkerSize',5);
        % COHERENCE
        plot(x,crespcoh_(ri,:),'-','Color',cmap(3,:));
%         errbar(coherences,cc_(ri,cohidx),squeeze(cc(2,ri,cohidx))-cc_(ri,cohidx)','-','Color',cmap(3,:));
        plot(coherences,squeeze(cc_(mi,ri,cohidx)),'o','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','white','MarkerSize',5);
        if any(ri==[1 2])
            axis([0 1 0 2.5]);
        else
            axis([0 1 0 1]);
        end
        if ri==7
            set(gca,'XTick',[0 1],'XTickLabeL',{'0%','100%'},'YTick',[0 1]);
            xlabel('Stimulus strength (%)');
            ylabel('\Delta signal (%)');
        else
            set(gca,'XTick',[0 1],'YTick',[0 1]);
        end
        title(sprintf('%s',ROIs{ri}));
        % TIMING
    %     subplot(8,2,(ri-1)*2+2); hold on
    %     plot(.250:.01:4,resptime_(ri,:),'-k');
    %     plot(times,time_(ri,timidx),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',3);
    %     
    %     axis([0.25 4 0 3]);
    %     if ri==8
    %         set(gca,'XTick',[0.5 1 2],'YTick',[0 1 2]);
    %         xlabel('Stimulus length (s)');
    %         ylabel('Signal change (%)');
    %     else
    %         set(gca,'XTick',[0.5 1 2],'YTick',[0 1 2]);
    %     end
        drawPublishAxis('figSize=[6,8.9]');
    end

    savepdf(h,fullfile(datafolder,'avg_models',sprintf('avg_behavtofmri_%s.pdf',models{mi})));
end