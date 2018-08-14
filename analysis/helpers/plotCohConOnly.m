function plotCohConOnly( nSIDs )
%PLOTCOHCONONLY Plot the cohcon deconvolution fit results
%
% This plots two things, a nice plot showing the average deconvolution and
% associated fits (should look reasonable)
%
% And separately it plots just the beta estimates across subjects for each
% condition, which should show 'ramping' behavior. 


%% Load Data

load(fullfile(datafolder,'avg_deconfits.mat'));
t = fits{1}.t;

rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

%% Move data into full datasets
adata.con = zeros(11,7,4,81);
adata.coh = zeros(11,7,5,81);
amodel.con = zeros(11,7,4,81);
amodel.coh = zeros(11,7,5,81);

stims = {'con','coh'};
for ni = 1:11
    for ri = 1:8
        adata.con(ni,ri,:,:) = squeeze(fits{ni}.data.con(ri,1,:,:));
        adata.coh(ni,ri,:,:) = squeeze(fits{ni}.data.coh(ri,1,:,:));
        amodel.con(ni,ri,:,:) = squeeze(fits{ni}.model.con(ri,1,:,:));
        amodel.coh(ni,ri,:,:) = squeeze(fits{ni}.model.coh(ri,1,:,:));
    end
end

%% Mean and std (do bootci later...)
adata.mcon = squeeze(mean(adata.con));
adata.mcoh = squeeze(mean(adata.coh));
amodel.mcon = squeeze(mean(amodel.con));
amodel.mcoh = squeeze(mean(amodel.coh));

%% Plot 1: Cohcon deconvolution plot
% See deconCohCon for similar plotting features
cmap = brewermap(11,'PuOr');

concolors = [-1 4 3 2];
cohcolors = [-1 7 8 9 10];
for ri = 1:8
    h = figure;
    subplot(211), hold on
    hs = zeros(1,3);
    for ci = 2:4
        hs(ci-1) = plot(t,squeeze(adata.mcon(ri,ci,:)),'o','MarkerFaceColor',cmap(concolors(ci),:),'MarkerEdgeColor',[1 1 1]);
        plot(t,squeeze(amodel.mcon(ri,ci,:)),'-','Color',cmap(concolors(ci),:));
    end
    legend(hs,{'+25%','+50%','+75%'});
    axis([0 30 -0.5 3])
    set(gca,'XTick',[0 10 20 30]);
    set(gca,'YTick',[-0.5 0 1]);
    xlabel('Time (s)');
    ylabel('Signal change (%)');
    drawPublishAxis
    subplot(212), hold on
    hs = zeros(1,4);
    for ci = 2:5
        hs(ci-1) = plot(t,squeeze(adata.mcoh(ri,ci,:)),'o','MarkerFaceColor',cmap(cohcolors(ci),:),'MarkerEdgeColor',[1 1 1]);
        plot(t,squeeze(amodel.mcoh(ri,ci,:)),'-','Color',cmap(cohcolors(ci),:));
    end
    legend(hs,{'+25%','+50%','+75%','+100%'});
    axis([0 30 -0.5 3])
    set(gca,'XTick',[0 10 20 30]);
    set(gca,'YTick',[-0.5 0 1]);
    xlabel('Time (s)');
    ylabel('Signal change (%)');
    drawPublishAxis
    savepdf(h,fullfile(datafolder,sprintf('avg_cc_nt_%s.pdf',rois{ri})));
end

%% Load and average fmri response data
% r2s = zeros(length(nSIDs),10,2);
% allcon = nan(length(nSIDs),101);
% allcoh = nan(length(nSIDs),101);
allconfmri = zeros(length(nSIDs),8,101);
allcohfmri = zeros(length(nSIDs),8,101);
offsets = zeros(length(nSIDs),8);
for ni = 1:length(nSIDs)
    sid = nSIDs(ni);
%     load(fullfile(datafolder,sprintf('s%03.0f_data.mat',sid)));
    load(fullfile(datafolder,sprintf('s%04.0f_fitroi.mat',sid)));
    x = 0:.01:1;
%     con = conModel(x,data.fits{1}.params);
%     coh = cohModel(x,data.fits{1}.params);
%     if ~isnan(sigmas.(sprintf('s%04.0f',sid))) && sigmas.(sprintf('s%04.0f',sid))<1 && sigmas.(sprintf('s%04.0f',sid))>0
%         allcon(ni,:) = con*sigmas.(sprintf('s%04.0f',sid));
%         allcoh(ni,:) = coh*sigmas.(sprintf('s%04.0f',sid));
%     end
    pos = 1:2:16;
    for ri = 1:8
        lroi = fitroi.roiparams{pos(ri)};
        rroi = fitroi.roiparams{pos(ri)+1};
        confmri = mean([conModel(x,lroi);conModel(x,rroi)]);
        cohfmri = mean([cohModel(x,lroi);cohModel(x,rroi)]);
        offsets(ni,ri) = mean([lroi.offset,rroi.offset]);
        allconfmri(ni,ri,:) = confmri-confmri(1);
        allcohfmri(ni,ri,:) = cohfmri-cohfmri(1);
%         r2s(ni,ri,1) = corr(coh',cohfmri');
%         r2s(ni,ri,2) = corr(con',confmri');
    end
end

fmri.con = allconfmri;
fmri.coh = allcohfmri;
%% Plot 2: Ramping plot

%
x = 0:.01:1;
ncon = squeeze(mean(fmri.con));
ncoh = squeeze(mean(fmri.coh));

rampcon = zeros(11,8,4);
rampcoh = zeros(11,8,5);
% move data into ramping datasets
for ni = 1:11
    rampcon(ni,:,:) = fits{ni}.roiparams.con;
    rampcoh(ni,:,:) = fits{ni}.roiparams.coh;
end
% mean data, and remove 0 condition
temp = bootci(1000,@mean,rampcon(:,:,2:end));
rampcon_ = squeeze(mean(temp));
rampcon_s = squeeze(temp(2,:,:,:))-rampcon_;
% rampcon_ = squeeze(mean(rampcon(:,:,2:end)));
temp = bootci(1000,@mean,rampcoh(:,:,2:end));
rampcoh_ = squeeze(mean(temp));
rampcoh_s = squeeze(temp(2,:,:,:))-rampcoh_;

con = [0.5 0.75 1];
coh = [0.25 0.5 0.75 1];
bcon = 0.25;
bcoh = 0;

moffsets = mean(offsets);

% interpolate the neural results
neuralcon_ = zeros(size(rampcon_));
neuralcoh_ = zeros(size(rampcoh_));
for ri = 1:length(rois)
    for ci = 1:length(con)
        tcon = con(ci)-bcon;
        neuralcon_(ri,ci) = 5*(ncon(ri,x==(tcon+bcon)) - ncon(ri,x==bcon)) + moffsets(ri);
    end
    for ci = 1:length(coh)
        tcoh = coh(ci)-bcoh;
        neuralcoh_(ri,ci) = 5*(ncoh(ri,x==(tcoh+bcoh)) - ncoh(ri,x==bcoh)) + moffsets(ri);
    end
end

h = figure;
cmap = brewermap(5,'PuOr');
subplot(211), hold on
plot(con,neuralcon_','-','Color',cmap(2,:));
plot(con,rampcon_','o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',[1 1 1]);
errbar(repmat(con,1,8),rampcon_',rampcon_s','Color',cmap(2,:));
set(gca,'XTick',[0.5 0.75 1],'XTickLabel',[50 75 100]);
set(gca,'YTick',[0 1 2 3]);

for ri = 1:8
    ypos = rampcon_(ri,end);
    if any(ri==[8 3])
        text(1.025,ypos,rois{ri});
    elseif any(ri==[4 6])
        text(1.045,ypos,rois{ri});
    else
        text(1.005,ypos,rois{ri});
    end
end

drawPublishAxis
subplot(212), hold on
plot(coh,neuralcoh_','-','Color',cmap(4,:));
plot(coh,rampcoh_','o','MarkerFaceColor',cmap(4,:),'MarkerEdgeColor',[1 1 1]);
errbar(repmat(coh,1,8),rampcoh_',rampcoh_s','Color',cmap(4,:));
xlabel('Stimulus strength (%)');
ylabel('Response amplitude (%)');
axis([0.25 1 0 1.25]);
set(gca,'XTick',[0.25 .5 1],'XTickLabel',[25 50 100]);
set(gca,'YTick',[0 1],'YTickLabeL',[0 1]);

for ri = 1:8
    ypos = rampcoh_(ri,end);
    if any(ri==[8 3 4])
        text(1.025,ypos,rois{ri});
    elseif any(ri==[4 6])
        text(1.035,ypos,rois{ri});
    else
        text(1.005,ypos,rois{ri});
    end
end
drawPublishAxis

savepdf(h,fullfile(datafolder,'avg_ramping.pdf'));
