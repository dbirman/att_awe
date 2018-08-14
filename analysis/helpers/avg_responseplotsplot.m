%% plot
load(fullfile(datafolder,'avg_fmriresponse.mat'));

f = figure; hold on
cmap = brewermap(7,'PuOr');
subplot(2,1,1), hold on
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
% xlabel('Contrast (%)');
% ylabel('Signal Change (%)');
set(gca,'XTick',[0 50 100]);
axis([0 100 0 0.8]);
set(gca,'XTick',[0 50 100]);
set(gca,'YTick',[0 0.2 0.4 0.6 0.8]);


subplot(2,1,2), hold on
% plot(x,mean(allcoh,1),'Color',cmap(6,:),'LineWidth',2);
plot(x*100,squeeze(mean(allcohfmri,1))','color',cmap(6,:));
axis([0 106 0 0.2]);
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
xlabel('Stimulus Strength (%)');
ylabel('Signal Change (%)');
set(gca,'XTick',[0 50 100]);
set(gca,'YTick',[0 0.2]);

savepdf(f,fullfile(datafolder,'avg_fmri','avg_fmriresponse.pdf'),'subplots=[2,1]','figsize=[6.5,10]');
