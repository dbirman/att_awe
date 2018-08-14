mcontrol = squeeze(nanmean(control,1));
controls = squeeze(nanstd(control,1));
mattend = squeeze(nanmean(attend,1));
attends = squeeze(nanstd(attend,1));
munattend = squeeze(nanmean(unattend,1));
unattends = squeeze(nanstd(unattend,1));

map = brewermap(6,'PuOr');
h = figure; hold on
h2 = plot([0.15 0.3 0.45 0.6],mcontrol(1,:),'-','Color',map(1,:));
    errbar([0.15 0.3 0.45 0.6],mcontrol(1,:),controls(1,:),'Color',map(1,:));
h6 = plot([0.325 0.4 0.55 0.85],mcontrol(2,:),'-','Color',map(6,:));
    errbar([0.325 0.4 0.55 0.85],mcontrol(2,:),controls(2,:),'Color',map(6,:));

h1 = plot([0.15 0.3 0.45 0.6],mcontrol(1,:),'o','MarkerSize',15);
set(h1(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(1,:),'LineWidth',1.5);

h3 = plot(0.3025,mattend(1),'o','MarkerSize',15);
set(h3(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(2,:),'LineWidth',1.5);
    errbar(0.3025,mattend(1),attends(1),'Color',map(2,:));
h7 = plot(0.305,munattend(1),'o','MarkerSize',15);
    errbar(0.305,munattend(1),unattends(1),'Color',map(3,:));
set(h7(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(3,:),'LineWidth',1.5);
h5 = plot([0.325 0.4 0.55 0.85],mcontrol(2,:),'o','MarkerSize',15);
set(h5(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(6,:),'LineWidth',1.5);
h9 = plot(0.4025,mattend(2),'o','MarkerSize',15);
    errbar(0.4025,mattend(2),attends(2),'Color',map(5,:));
set(h9(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(5,:),'LineWidth',1.5);

h4 = plot(0.405,munattend(2),'o','MarkerSize',15);
    errbar(0.405,munattend(2),unattends(2),'Color',map(4,:));
set(h4(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(4,:),'LineWidth',1.5);
%    ;


legend([h1,h3,h7,h5,h9,h4],{'Coherence: Control','Coherence: Attended','Coherence: Unattended','Contrast: Control','Contrast: Attended','Contrast: Unattended'});

title('Psychometric Functions for Cohcon');
xlabel('Contrast/Coherence (%)');
ylabel('Threshold (%)');
axis([0.1 0.9 0 1])
drawPublishAxis

fname = fullfile('C:/Users/Dan/Documents/Box Sync/COHCON_DATA','avg_performance.pdf');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'InvertHardCopy','off');
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fname,'-dpdf');