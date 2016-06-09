function subj_behav_analysis( subj, refit )
%%
cfolder = sprintf('~/data/cohcon/%s',subj);
% cfolder = '~/data/cohcon/s300';

ffolder = fullfile(cfolder,'figures');
if ~isdir(ffolder), mkdir(ffolder); end


files = dir(fullfile(cfolder,'*.mat'));
%% dispInfo
fname = fullfile(ffolder,'performance.pdf');
if ~isfile(fname)
    load(fullfile(cfolder,files(end).name));
    h = ccDispInfo(stimulus,subj);
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'InvertHardCopy','off');
    set(gcf,'Color',[1 1 1]);
    set(gca,'Color',[1 1 1]);
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fname,'-dpdf');
end

%% Load data

% format:
%     1       2         3        4      5      6     7      8       9
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch -
%      10      11        12
%   pedcon - pedcoh - correct

adata = [];
for fi = 1:length(files)
    load(fullfile(cfolder,files(fi).name));
    e = getTaskParameters(myscreen,task);
    if length(e{1})==2
        e = e{1}(2);
        flip = [-1 1];
        adata = [adata ; [repmat(stimulus.runs.curTask*flip(stimulus.nocatch+1),length(e.response),1), ...
            repmat(stimulus.baseCon,length(e.response),1), repmat(stimulus.baseCoh,length(e.response),1), ...
            e.randVars.lCon', e.randVars.rCon',...
            e.randVars.lCoh', e.randVars.rCoh',...
            (e.response-1)', (e.parameter.catch)',...
            e.randVars.contrast', e.randVars.coherence',e.randVars.correct']];
    end
end

%% Remove NaN
adata = adata(~any(isnan(adata),2),:);

%% Fit behav model
if refit
    fit = fitCCBehavModel(adata,1,'con-naka,coh-linear');

    %% Try other models
    fits.allnaka = fitCCBehavModel(adata,1,'con-naka,coh-naka');
    fits.alllinear = fitCCBehavModel(adata,1,'con-linear,coh-linear');
    % Options:
    % noalpha - doesn't fit the mixture model for split attention
    % nounatt - doesn't reduce response intensity when unattended
    % unattnoise - fits the unattended condition using larger noise
    % poisson - poisson noise rather than regular noise
    fits.nobias = fitCCBehavModel(adata,1,'con-naka,coh-linear,nobias');
    fits.nounatt = fitCCBehavModel(adata,1,'con-naka,coh-linear,nounatt');
    % fits.unattnoise = fitCCBehavModel(adata,1,'con-naka,coh-linear,unattnoise');
    fits.poisson = fitCCBehavModel(adata,1,'con-naka,coh-linear,poisson');
end

%% Save data
load(sprintf('~/data/cohcon/%s_data.mat',subj));
if refit
    data.fit = fit;
    data.fits = fits;
else
    fit = data.fit;
end
save(sprintf('~/data/cohcon/%s_data.mat',subj),'data');

%% Figure
x = 0:.01:1;
con = x;
cony = fit.params.conRmax .* ((con.^fit.params.conn) ./ (con.^fit.params.conn + fit.params.conc50.^fit.params.conn));
cony_un = (fit.params.conRmax * fit.params.conunatt) .* ((con.^fit.params.conn) ./ (con.^fit.params.conn + fit.params.conc50.^fit.params.conn));

cohy = fit.params.cohslope*x;
cohy_un = fit.params.cohslope*x*fit.params.cohunatt;

h=figure;
hold on
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
fname = fullfile(ffolder,'response.pdf');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'InvertHardCopy','off');
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fname,'-dpdf');

%%
h = cc_rightchoice(adata, fit);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
title(sprintf('BIC: %4.2f',fit.BIC));
fname = fullfile(ffolder,'rightchoice.pdf');
set(h,'InvertHardCopy','off');
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
print(fname,'-dpdf');
%% Thresholds over time
% thresholdTimePlot(stimulus);

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