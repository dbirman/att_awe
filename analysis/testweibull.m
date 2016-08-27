function testweibull(subj)

%% Testing how different weibull fits are:
if strfind(ver('os'),'Windows')
    load(sprintf('/Users/dan/proj/Box Sync/COHCON_DATA/%s_data.mat',subj));
else
    load(fullfile(sprintf('~/Desktop/Box Sync/COHCON_DATA/%s_data.mat',subj)));
end
alpha = data.fit.params.conalpha;
kappa = data.fit.params.conkappa;

x = 0:.001:2;
y = -alpha*exp(-kappa*x)+alpha;

%%
pedestals = [0.325 0.4 0.55 0.85];
truetruethresholds = zeros(size(pedestals));
for pi = 1:length(pedestals)
    pedestalsR(pi) = y(x==pedestals(pi));
    truethresholds(pi) = 1./(alpha*kappa*exp(-kappa*pedestals(pi)));
    cy = pedestalsR(pi);
    truetruethresholds(pi) = x(find(y>=(cy+1),1))-pedestals(pi);
end

%% simulate with sigma = 1

thresholds = zeros(size(pedestals));
thresholds_ = zeros(size(pedestals));

numR = 1000;
range = [0:.005:0.15 0.2 0.3];
for pi = 1:length(pedestals)
    thisrange = zeros(size(range));
    for ri = 1:length(range)
        correct = zeros(1,numR);
        for i = 1:numR
            if (randn+pedestalsR(pi)) < (randn+y(find(x>=(pedestals(pi)+range(ri)),1)))
                correct(i) = 1;
            end
        end
        thisrange(ri) = sum(correct);
    end
    out = fitweibull(range,thisrange,'ntotal',repmat(numR,1,length(range)));
    thresholds(pi) = out.fitparams(1);
    thresholds_(pi) = out.x(find(out.y>=0.76,1));
end

%% Plot
f = figure; hold on

cmap = brewermap(6,'PuOr');
h(1) = plot(pedestals,data.control(2,:),'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',15);
h(2) = plot(pedestals,thresholds,'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',15);
h(3) = plot(pedestals+.01,thresholds_,'-o','Color',cmap(3,:),'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',15);
h(4) = plot(pedestals,truethresholds,'-o','Color',cmap(4,:),'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',15);
h(5) = plot(pedestals,truetruethresholds,'-o','Color',cmap(5,:),'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',15);
legend(h,{'Actual Thresholds','Simulated @ 82%','Simulated @ 76%','+1 Response Estimate','1/Derivative Estimate'});
xlabel('Contrast Pedestal (%)');
ylabel('Just Noticeable Difference (%)');
drawPublishAxis
h = figure(f);
fname = fullfile('~/Desktop/hella_confused.pdf');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'InvertHardCopy','off');
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fname,'-dpdf');
