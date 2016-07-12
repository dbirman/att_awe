function subj_behav_analysis( subj, refit )

disp(sprintf('Started running for %s',subj));
%%

%% Load data
adata = loadadata(subj);

%% Fit Contrast/Coherence response models (just to control condition)
out = fitCCBehavControlModel(adata);
%% Fit behav model
fit = struct;
fits = struct;
if refit
    fit = fitCCBehavModel(adata,0,'con-naka,coh-linear,nounatt'); % uses beta weights to fit the unattended conditions

    %% Try other models
    fits.allnaka = fitCCBehavModel(adata,0,'con-naka,coh-naka,nounatt');
    fits.alllinear = fitCCBehavModel(adata,0,'con-linear,coh-linear,nounatt');
    fits.nobias = fitCCBehavModel(adata,0,'con-naka,coh-linear,nounatt,nobias');
    fits.unattnoise = fitCCBehavModel(adata,0,'con-naka,coh-linear,unattnoise');
    fits.poisson = fitCCBehavModel(adata,0,'con-naka,coh-linear,nounatt,poisson');
end

%% Save data
if exist(fullfile(sprintf('C:/Users/Dan/proj/COHCON_DATA/%s_data.mat',subj)))==2
    load(fullfile(sprintf('C:/Users/Dan/proj/COHCON_DATA/%s_data.mat',subj)));
end
if refit
    data.fit = fit;
    data.fits = fits;
else
    fit = data.fit;
end
save(fullfile(sprintf('C:/Users/Dan/proj/COHCON_DATA/%s_data.mat',subj)),'data');

return

%% dispInfo
load(fullfile(cfolder,files(end).name));
h1 = ccDispInfo(stimulus,subj);
fname = fullfile('C:/Users/Dan/proj/COHCON_DATA/',sprintf('%s_threshold.pdf',subj));

set(h1,'Units','Inches');
pos = get(h1,'Position');
set(h1,'InvertHardCopy','off');
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fname,'-dpdf');
%% Figure
h1 = figure;
x = 0:.01:1;
con = x;
cony = fit.params.conRmax .* ((con.^fit.params.conn) ./ (con.^fit.params.conn + fit.params.conc50.^fit.params.conn));
cony_un = (fit.params.conRmax * fit.params.conunatt) .* ((con.^fit.params.conn) ./ (con.^fit.params.conn + fit.params.conc50.^fit.params.conn));

cohy = fit.params.cohslope*x;
cohy_un = fit.params.cohslope*x*fit.params.cohunatt;

clist = brewermap(3,'PuOr');
plot(x,cony,'Color',clist(1,:));
plot(x,cohy,'Color',clist(3,:));
if fit.params.poissonNoise
    [hl,hp] = boundedline(x,cony,sqrt(fit.params.sigma*cony),'alpha');
else
    [hl,hp] = boundedline(x,cony,fit.params.sigma*ones(size(cony)),'alpha');
end
set(hl,'Color',clist(1,:));
set(hp,'FaceColor',clist(1,:));
[hl,hp] = boundedline(x,cony_un,fit.params.sigma*ones(size(cony)),'alpha');
set(hl,'Color',clist(1,:));
set(hp,'FaceColor',clist(1,:));
[hl,hp] = boundedline(x,cohy,fit.params.sigma*ones(size(cohy)),'alpha');
set(hl,'Color',clist(3,:));
set(hp,'FaceColor',clist(3,:));
[hl,hp] = boundedline(x,cohy_un,fit.params.sigma*ones(size(cohy)),'alpha');
set(hl,'Color',clist(3,:));
set(hp,'FaceColor',clist(3,:));
legend({'Contrast','Coherence'});
xlabel('Contrast / Coherence (%)');
ylabel('Response (a.u.)');

drawPublishAxis
% fname = fullfile(ffolder,'response.pdf');
% set(h2,'Units','Inches');
% pos = get(h2,'Position');
% set(h2,'InvertHardCopy','off');
% set(gcf,'Color',[1 1 1]);
% set(gca,'Color',[1 1 1]);
% set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fname,'-dpdf');
fname = fullfile('C:/Users/Dan/proj/COHCON_DATA/',sprintf('%s_response.pdf',subj));

set(h1,'Units','Inches');
pos = get(h1,'Position');
set(h1,'InvertHardCopy','off');
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fname,'-dpdf');


%%
h1 = cc_rightchoice(adata, fit);
% set(h3,'Units','Inches');
% pos = get(h3,'Position');
% set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% title(sprintf('BIC: %4.2f',fit.BIC));
% fname = fullfile(ffolder,'rightchoice.pdf');
% set(h3,'InvertHardCopy','off');
% set(gcf,'Color',[1 1 1]);
% set(gca,'Color',[1 1 1]);
% print(fname,'-dpdf');

%%
fname = fullfile('C:/Users/Dan/proj/COHCON_DATA/',sprintf('%s_rightchoice.pdf',subj));

set(h1,'Units','Inches');
pos = get(h1,'Position');
set(h1,'InvertHardCopy','off');
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fname,'-dpdf');

%%
disp(sprintf('Finished running for %s',subj));

%% Check for stay/switch bias
% generate a matrix that is:
%   P(right | prevR/prevL, succ/fail)
%
%         prevFail  prevSucc
%   prevR    P        P
%   prevL    P        P


% format:
%     1       2         3        4      5      6     7      8       9
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch -
%      10      11        12
%   pedcon - pedcoh - correct
% mat = zeros(2,2);
% n = zeros(2,2);
% for succ = 0:1 % fail/succ
%     for prev = 0:1 % left/right
%         for ai = 2:size(adata,1)
%             pdat = adata(ai-1,:);
%             dat = adata(ai,:);
%             if pdat(12) == succ && pdat(8)==prev
%                 % use this trial
%                 mat(prev+1,succ+1) = mat(prev+1,succ+1)+ dat(8);
%                 n(prev+1,succ+1) = n(prev+1,succ+1)+1;
%             end
%         end
%     end
% end
% mat./n