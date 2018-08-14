function [ycon,ycoh] = avg_behavmodel( sids, fitnum )

x = 0:.01:1;
ycon = zeros(length(sids),length(x));
ycoh = zeros(length(sids),length(x));

for si = 1:length(sids)
    load(fullfile(datafolder,sprintf('%s_data.mat',sids{si})));
    ycon(si,:) = conModel(x,data.fits{fitnum}.params);
    ycoh(si,:) = cohModel(x,data.fits{fitnum}.params);
end

h = figure; hold on
cmap = brewermap(7,'PuOr');
for i = 1:length(sids)
    plot(x,ycon(i,:),'-','Color',cmap(2,:),'LineWidth',0.25);
    plot(x,ycoh(i,:),'-','Color',cmap(6,:),'LineWidth',0.25);
end
plot(x,mean(ycon,1),'-','Color',cmap(2,:),'LineWidth',3);
plot(x,mean(ycoh,1),'-','Color',cmap(6,:),'LineWidth',3);
xlabel('Stimulus strength (%)');
set(gca,'XTick',[0 0.5 1],'XTickLabel',[0 50 100]);
ylabel('Model response (a.u.)');
set(gca,'YTick',[0 30],'YTickLabel',{'',''});
if mod(fitnum,2)==1
    axis([0 1 0 30]);
else
    axis([0 1 0 140]);
end

p(1) = plot(-1,-1,'-k','LineWidth',3);
p(2) = plot(-1,-1,'-k','LineWidth',1);
legend(p,{'Average','Individual'});

strs = {'con-exp,coh-exp','con-exp,coh-exp,poisson','con-linear,coh-linear','con-linear,coh-linear,poisson','con-naka,coh-naka','con-naka,coh-naka,poisson','con-explin,coh-explin','con-explin,coh-explin,poisson'}; %
type = strs{fitnum};
drawPublishAxis('figSize=[8.5,10]');

savepdf(h,fullfile(datafolder,'avg_behav',sprintf('avg_behavmodel_%s.pdf',type)));

%% Plot the model figure while we're here

% h = figure; hold on
% mx = 0:.01:1;
% my = mean(ycon,1);
% plot(x,my,'-k','LineWidth',3);
% xlabel('Stimulus strength (%)');
% set(gca,'XTick',[0 0.5 1],'XTickLabel',[0 50 100]);
% ylabel('Model response (s.d.)');
% set(gca,'YTick',[0 30],'YTickLabel',{'',''});
% 
% vals = [15 70];
% for vi = 1:length(vals)
%     plot([mx(vals(vi)) mx(vals(vi))],[0 my(vals(vi))],'--k');
%     yval = my(vals(vi));
%     plot([0 mx(vals(vi))],[yval yval],'--k');
%     yval2 = yval+2.5;
%     xval2 = find(my>=yval2,1);
%     plot([0 mx(xval2)],[yval2 yval2],'--k');
%     plot([mx(xval2) mx(xval2)],[0 yval2],'--k');
%     text(-.075,yval+(yval2-yval)/2,'+1 s.d.');
%     text(mx(vals(vi))+(mx(xval2)-mx(vals(vi)))/2-0.055,-0.75,'JND');
% end
% 
% savepdf(h,fullfile(datafolder,'avg_behav','display_model.pdf'),'figsize=[8,10]');

