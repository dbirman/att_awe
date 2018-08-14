%%
nSIDs = [   305   329    43    25   300   346   343   344   338   340   348];


%% Fit different exponent models

expfits = cell(size(nSIDs));
hrfexp = zeros(length(nSIDs),8);
for si = 1:length(nSIDs)
    decondata = decondatas{si}.decondata;
    allData = reformAllData(decondata);
    
    load(fullfile(datafolder,sprintf('s%04.0f_hrf.mat',nSIDs(si)))); % saved as cur
    expfits{si} = fitCCHRFModel_exp(allData,'fitexp',cur);
    for ri = 1:8
        hrfexp(si,ri) = expfits{si}.roifit{ri}.params.hrfexp;
    end
    close all
end

%% Plot HRF exp across ROIs
hrfexp_ = squeeze(bootci(1000,@mean,hrfexp));
hrfexp = squeeze(mean(hrfexp_));

%% Plot
h = figure; hold on
errbar(1:8,hrfexp,squeeze(hrfexp_(2,:))-hrfexp,'k');
plot(1:8,hrfexp,'o','MarkerFaceColor','k','MarkerEdgeColor','w');
set(gca,'XTick',1:8,'XTickLabel',rois);
ylabel('Decay exponent');
axis([1 8 -0.675 -0.575]);
set(gca,'YTick',[-0.65 -0.6]);
drawPublishAxis('figSize=[8.5 4.25]');
savepdf(h,fullfile(datafolder,'avg_fmri','expfit.pdf'));

%% BIC comparison
figure
for si = 1:length(nSIDs)
    subplot(length(nSIDs),1,si); hold on
%     plot(fits{si,1}.BIC,'g');
    hline(0,'g');
    plot(fits{si,2}.BIC-fits{si,1}.BIC,'y');
    plot(fits{si,3}.BIC-fits{si,1}.BIC,'r');
    drawPublishAxis
end

%% Load
% subject | roi | model | index | value
time = nan(length(nSIDs),length(ROIs),2,20,81);
cc = nan(length(nSIDs),length(ROIs),2,20,81);
for ni = 1:length(nSIDs)
    for ri = 1:length(ROIs)
        if ~any(ni==[3 4 5])
            time(ni,ri,1,:,:) = fits{ni,1}.roifit{ri}.time.resp;
            time(ni,ri,2,:,:) = fits{ni,1}.roifit{ri}.time.model;
        end
        cc(ni,ri,1,:,:) = fits{ni,1}.roifit{ri}.cc.resp;
        cc(ni,ri,2,:,:) = fits{ni,1}.roifit{ri}.cc.model;
    end
end

%% Average over subjects
time_ = nan(2,2,20,81);
time_s = time_;
cc_ = time_;
cc_s = time_;
reps = 250;
disppercent(-1/length(rois));
for ri = 1:length(rois)
    for mi = 1:2
        for ii = 1:20
            tout = bootci(reps,@nanmean,squeeze(time(:,ri,mi,ii,:)));
            toutm = squeeze(mean(tout));
            time_(ri,mi,ii,:) = toutm;
            time_s(ri,mi,ii,:) = tout(2,:)-toutm;
            cout = bootci(reps,@nanmean,squeeze(cc(:,ri,mi,ii,:)));
            coutm = squeeze(mean(cout));
            cc_(ri,mi,ii,:) = coutm;
            cc_s(ri,mi,ii,:) = cout(2,:)-coutm;
        end
    end
    disppercent(ri/length(rois));
end
disppercent(inf);

%%
avg_decon.time = time_;
avg_decon.times = time_s;
avg_decon.cc = cc_;
avg_decon.ccs = cc_s;
avg_decon.conidx = decondata.V1.cc.conidxs;
avg_decon.cohidx = decondata.V1.cc.cohidxs;
avg_decon.timconidx = decondata.V1.time.conidxs;
avg_decon.timcohidx = decondata.V1.time.cohidxs;
avg_decon.timidxs = decondata.V1.time.timidxs;

save(fullfile(datafolder,'avg_decon_deconv.mat'),'avg_decon');

%% convert fits to con/coh models
avg_conv = struct;
avg_conv.con = zeros(11,8,101); avg_conv.coh = zeros(11,8,101);

for si = 1:length(nSIDs)
    cfit = fits{si,1};
    
    x = 0:.01:1;
    for ri = 1:8
        avg_conv.con(si,ri,:) = conModel(x,cfit.roifit{ri}.params);
        avg_conv.coh(si,ri,:) = cohModel(x,cfit.roifit{ri}.params);
    end
end


%% plot allconfmri/allcohfmri
x = 0:.01:1;
h = figure; hold on
cmap = brewermap(7,'PuOr');
for si = 1:length(nSIDs)
    plot(x,squeeze(avg_conv.con(si,1,:)),'-','Color',cmap(2,:));
    plot(x,squeeze(avg_conv.coh(si,8,:)),'-','Color',cmap(6,:));
end
conv1 = mean(bootci(1000,@nanmean,squeeze(avg_conv.con(:,1,:))));
cohmt = mean(bootci(1000,@nanmean,squeeze(avg_conv.coh(:,8,:))));
plot(x,conv1,'Color',cmap(2,:),'LineWidth',4);
plot(x,cohmt,'Color',cmap(6,:),'LineWidth',4);
axis([0 1 0 2.5])
set(gca,'XTick',[0 .5 1],'XTickLabel',[0 50 100]);
set(gca,'YTick',[0 1 2]);
xlabel('Stimulus Strength (%)');
ylabel('Signal Change (%)');
drawPublishAxis
savepdf(h,fullfile(datafolder,'avg_fmrifunctions_conv.pdf'));

%% compute average offset parameter

offsets = zeros(11,8);
for si = 1:length(nSIDs)
    for ri = 1:length(ROIs)
        offsets(si,ri) = fits{si,1}.roifit{ri}.params.offset;
    end
end
ci = bootci(1000,@nanmean,offsets);
mu = mean(ci,1);

h = figure; hold on

plot(1:8,mu,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1],'MarkerSize',15);
errbar(1:8,mu,ci(2,:)-mu,'-k');

set(gca,'XTick',1:8,'XTickLabel',ROIs);
xlabel('ROI');
ylabel('Offset effect (% signal change)');

drawPublishAxis

savepdf(h,fullfile(datafolder,'avg_offset.pdf'));