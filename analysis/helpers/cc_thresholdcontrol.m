function cc_thresholdcontrol( subj )
%UNTITLED Estimate the threshold values for each output (depends only on
%the contrast/coherence response functions, doesn't need to deal with the
%beta values at all). Basically we're inverting the signal detection model
%directly.


%     1       2         3        4      5      6     7      8       9
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch -
%      10      11        12
%   pedcon - pedcoh - correct

if strfind(getenv('OS'),'Windows')
    load(fullfile(sprintf('C:/Users/Dan/proj/COHCON_DATA/%s_data.mat',subj)));
else
    warning('not implemented');
end

adata = loadadata(subj);

fit = data.fit;

%%
x = 0:.001:3;

% compute for the initial fit model
con = conModel(x,fit.params);
coh = cohModel(x,fit.params);

%% figure out the con/coh pedestals
conped = unique(adata(:,10));
cohped = unique(adata(:,11));

%% compute the thresholds for contrast
cont = zeros(size(conped));
for ci = 1:length(conped)
    cped = conped(ci);
    baseResp = conModel(cped,fit.params);
    % threshold = sigma
    if fit.params.poissonNoise
        cont(ci) = sqrt(fit.params.sigma*cped);
    else
        cont(ci) = fit.params.sigma;
    end
    % go find this value 
    cont(ci) = x(find(con>=(cont(ci)+baseResp),1))-cped;
end

%% compute for coherence
coht = zeros(size(cohped));
for ci = 1:length(cohped)
    cped = cohped(ci);
    baseResp = cohModel(cped,fit.params);
    if fit.params.poissonNoise
        noise = sqrt(fit.params.sigma*cped);
    else
        noise = fit.params.sigma;
    end
    coht(ci) = x(find(coh>=noise+baseResp,1))-cped;
end

data.control(data.control<=0) = NaN;
%% plot

map = brewermap(6,'PuOr');
h = figure; hold on

plot(conped,cont,'-','Color',map(1,:));
plot(cohped,coht,'-','Color',map(6,:));
h1 = plot(conped,data.control(2,:),'o','MarkerSize',15);
set(h1(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(1,:),'LineWidth',1.5);
h1 = plot(cohped,data.control(1,:),'o','MarkerSize',15);
set(h1(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(6,:),'LineWidth',1.5);

title(sprintf('Model Fit for Naka-Con, Linear-Coh, no Bias, Subject: %s',subj));
xlabel('Contrast/Coherence Pedestal (%)');
ylabel('Just Noticeable Difference (% Con/Coh)');
drawPublishAxis


fname = fullfile('C:/Users/Dan/proj/COHCON_DATA',sprintf('%s_thresholdmodel.pdf',subj));
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'InvertHardCopy','off');
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fname,'-dpdf');
%%
data.cont = cont;
data.coht = coht;

if strfind(getenv('OS'),'Windows')
    save(fullfile(sprintf('C:/Users/Dan/proj/COHCON_DATA/%s_data.mat',subj)),'data');
else
    warning('not implemented');
end