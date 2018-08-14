
% load(fullfile(datafolder,'sigmas.mat'));
r2s = zeros(length(nSIDs),10,2);
allcon = nan(length(nSIDs),101);
allcoh = nan(length(nSIDs),101);
allconfmri = zeros(length(nSIDs),8,101);
allcohfmri = zeros(length(nSIDs),8,101);
for ni = 1:length(nSIDs)
    sid = nSIDs(ni);
    load(fullfile(datafolder,sprintf('s%03.0f_data.mat',sid)));
    load(fullfile(datafolder,sprintf('s%04.0f_fitroi.mat',sid)));
    x = 0:.01:1;
    con = conModel(x,data.fits{1}.params);
    coh = cohModel(x,data.fits{1}.params);
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
        allconfmri(ni,ri,:) = confmri-confmri(1);
        allcohfmri(ni,ri,:) = cohfmri-cohfmri(1);
%         r2s(ni,ri,1) = corr(coh',cohfmri');
%         r2s(ni,ri,2) = corr(con',confmri');
    end
end
%% plot

warning('this code doesn''t really make sense... you''re fitting functions to each other, of course its going to look good.');
f = figure;
cmap = brewermap(7,'PuOr');
subplot(211), hold on
% plot(x,mean(allcon,1),'Color',cmap(2,:),'LineWidth',2);
plot(x*100,squeeze(mean(allconfmri,1))','Color',cmap(2,:));
axis([0 100 0 0.6]);
% annotate with rois
maf = squeeze(mean(allconfmri,1));
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
for ri = 1:8
    ypos = maf(ri,end);
    if any(ri==[3 4 6 8])
        text(106,ypos,rois{ri});
    else
        text(101,ypos,rois{ri});
    end
end
xlabel('Contrast (%)');
ylabel('Signal Change (%)');
drawPublishAxis
subplot(212), hold on
% plot(x,mean(allcoh,1),'Color',cmap(6,:),'LineWidth',2);
plot(x*100,squeeze(mean(allcohfmri,1))','color',cmap(6,:));
axis([0 100 0 0.15]);
% annotate with rois
maf = squeeze(mean(allcohfmri,1));
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
for ri = 1:8
    ypos = maf(ri,end);
    if ri==6
        text(106,ypos,rois{ri});
    else
        text(101,ypos,rois{ri});
    end
end
xlabel('Motion Coherence (%)');
ylabel('Signal Change (%)');
drawPublishAxis
savepdf(f,fullfile(datafolder,'avg_fmriresponse.pdf'));