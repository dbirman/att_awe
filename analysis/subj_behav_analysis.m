function subj_behav_analysis( subj, modes )

disp(sprintf('Started running for %s',subj));
%% Get files
files = dir(fullfile(datafolder,subj));

%% Load data
adata = loadadata(subj);

%% TEMP CODE
%% load original models
load(fullfile(datafolder,sprintf('%s_data.mat',subj)));
fits = data.fits;
BICs = data.BICs;

%% Fit Models
if strfind(modes,'refit')
    %% Fit Contrast/Coherence response models (just to control condition)
    strs = {'con-exp,coh-exp','con-exp,coh-exp,poisson','con-linear,coh-linear','con-linear,coh-linear,poisson'}; %
    %fits = cell(1,length(strs));
%     BICs = zeros(size(fits));
    minl = inf;
    for si = 3:4
        fits{si} = fitCCBehavControlModel(adata,1,strs{si});
        BICs(si) = fits{si}.BIC;
        if fits{si}.BIC < (minl-5)
            minl = fits{si}.BIC;
            fit = fits{si};
        end
    end
end
%% Save data
if exist(fullfile(datafolder,sprintf('%s_data.mat',subj)))==2
    load(fullfile(datafolder,sprintf('%s_data.mat',subj)));
end
if strfind(modes,'refit')
    data.fit = fit;
    data.fits = fits;
    data.BICs = BICs;
    data.strs = strs;
else
    fit = data.fit;
end
save(fullfile(datafolder,sprintf('%s_data.mat',subj)),'data');
% 

%%
if strfind(modes,'disp')
    %% dispInfo
    load(fullfile(datafolder,sprintf('%s',subj),files(end).name));
    ccDispInfocontrol(stimulus,subj);
%     fname = fullfile(datafolder,sprintf('%s_threshold.pdf',subj));

%     set(h1,'Units','Inches');
%     pos = get(h1,'Position');
%     set(h1,'InvertHardCopy','off');
%     set(gcf,'Color',[1 1 1]);
%     set(gca,'Color',[1 1 1]);
%     set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(fname,'-dpdf');

    %% Figure
%     h1 = figure;
%     x = 0:.01:1;
%     con = x;
%     cony = conModel(x,fit.params);
%     cony_un = (fit.params.conRmax * fit.params.conunatt) .* ((con.^fit.params.conn) ./ (con.^fit.params.conn + fit.params.conc50.^fit.params.conn));
% 
%     cohy = cohModel(x,fit.params);
%     cohy = fit.params.cohslope*x;
%     cohy_un = fit.params.cohslope*x*fit.params.cohunatt;
% 
%     clist = brewermap(3,'PuOr');
%     plot(x,cony,'Color',clist(1,:));
%     plot(x,cohy,'Color',clist(3,:));
%     if fit.params.poissonNoise
%         [hl,hp] = boundedline(x,cony,sqrt(fit.params.sigma*cony),'alpha');
%     else
%         [hl,hp] = boundedline(x,cony,fit.params.sigma*ones(size(cony)),'alpha');
%     end
%     set(hl,'Color',clist(1,:));
%     set(hp,'FaceColor',clist(1,:));
%     [hl,hp] = boundedline(x,cony_un,fit.params.sigma*ones(size(cony)),'alpha');
%     set(hl,'Color',clist(1,:));
%     set(hp,'FaceColor',clist(1,:));
%     [hl,hp] = boundedline(x,cohy,fit.params.sigma*ones(size(cohy)),'alpha');
%     set(hl,'Color',clist(3,:));
%     set(hp,'FaceColor',clist(3,:));
%     [hl,hp] = boundedline(x,cohy_un,fit.params.sigma*ones(size(cohy)),'alpha');
%     set(hl,'Color',clist(3,:));
%     set(hp,'FaceColor',clist(3,:));
%     legend({'Contrast','Coherence'});
%     xlabel('Contrast / Coherence (%)');
%     ylabel('Response (a.u.)');
% 
%     drawPublishAxis
%     % fname = fullfile(ffolder,'response.pdf');
%     % set(h2,'Units','Inches');
%     % pos = get(h2,'Position');
%     % set(h2,'InvertHardCopy','off');
%     % set(gcf,'Color',[1 1 1]);
%     % set(gca,'Color',[1 1 1]);
%     % set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     % print(fname,'-dpdf');
%     fname = fullfile(datafolder,sprintf('%s_response.pdf',subj));
% 
%     set(h1,'Units','Inches');
%     pos = get(h1,'Position');
%     set(h1,'InvertHardCopy','off');
%     set(gcf,'Color',[1 1 1]);
%     set(gca,'Color',[1 1 1]);
%     set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(fname,'-dpdf');
end
%%
if strfind(modes,'right')
    h1 = cc_rightchoicecontrol(adata, fit, 0);
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
    savepdf(h1,fullfile(datafolder,sprintf('%s_rightchoice.pdf',subj)));
end
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