function avg_comparebehavmodels( sids )
%%
x = 0:.01:1;
% exponentials
e.ycon = zeros(length(sids),length(x));
e.ycoh = zeros(length(sids),length(x));
% poisson
p.ycon = zeros(length(sids),length(x));
p.ycoh = zeros(length(sids),length(x));

for si = 1:length(sids)
    load(fullfile(datafolder,sprintf('%s_data.mat',sids{si})));
    e.ycon(si,:) = conModel(x,data.fits{1}.params);
    e.ycoh(si,:) = cohModel(x,data.fits{1}.params);
    p.ycon(si,:) = conModel(x,data.fits{7}.params);
    p.ycoh(si,:) = cohModel(x,data.fits{7}.params);
end

h = figure; hold on
cmap = brewermap(7,'PuOr');

plot(-1,-1,'-k');
plot(-1,-1,'--k');
legend({'Nonlinear additive','Linear poisson'});
boundedline(x,mean(p.ycon,1),sqrt(mean(p.ycon,1)),'--','cmap',cmap(2,:));
plot(x,mean(p.ycon,1),'--','Color',cmap(2,:),'LineWidth',1);
boundedline(x,mean(p.ycoh,1),sqrt(mean(p.ycoh,1)),'--','cmap',cmap(6,:));
plot(x,mean(p.ycoh,1),'--','Color',cmap(6,:),'LineWidth',1);
boundedline(x,mean(e.ycon,1),1,'cmap',cmap(2,:));
plot(x,mean(e.ycon,1),'-','Color',cmap(2,:),'LineWidth',1);
boundedline(x,mean(e.ycoh,1),1,'cmap',cmap(6,:));
plot(x,mean(e.ycoh,1),'-','Color',cmap(6,:),'LineWidth',1);
xlabel('\Delta stimulus (%)');
set(gca,'XTick',[0 0.5 1],'XTickLabel',[0 50 100]);
ylabel('Model response (s.d.)');
set(gca,'YTick',[0 20],'YTickLabel',{'',''});
axis([0 1 0 140]);
drawPublishAxis('figSize=[8.5,8]');
% types = {'Additive','Poisson'};
savepdf(h,fullfile(datafolder,'avg_behav','avg_behavmodel_compare.pdf'));

