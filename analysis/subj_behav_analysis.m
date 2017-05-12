function subj_behav_analysis( subj, modes )

disp(sprintf('Started running for %s',subj));
%% Get files
files = dir(fullfile(datafolder,subj));

%% Load data
adata = loadadata(subj);

%% Compute lapse rate
tdata = adata(adata(:,9)==-1,:);
tdata(:,13) = abs(tdata(:,5)-tdata(:,4)); % contrast diff
tdata(:,14) = abs(tdata(:,7)-tdata(:,6)); % coherence diff
conbins = quantile(tdata(:,13),.95);
cohbins = quantile(tdata(:,14),.95);
conidxs = logical((tdata(:,1)==2) .* (tdata(:,13)>conbins));
cohidxs = logical((tdata(:,1)==1) .* (tdata(:,14)>cohbins));
lapse = 1-mean([tdata(conidxs,12) ; tdata(cohidxs,12)]);

%% TEMP CODE
%% load original models
load(fullfile(datafolder,sprintf('%s_data.mat',subj)));
fits = data.fits;
BICs = data.BICs;

%% Fit Models
if strfind(modes,'refit')
    %% Fit Contrast/Coherence response models (just to control condition)
    strs = {'con-exp,coh-exp','con-exp,coh-exp,poisson','con-linear,coh-linear','con-linear,coh-linear,poisson','con-naka,coh-naka','con-naka,coh-naka,poisson','con-explin,coh-explin','con-explin,coh-explin,poisson'}; %
    %fits = cell(1,length(strs));
%     BICs = zeros(size(fits));
    minl = inf;
    for si = 1:length(strs)
        fits{si} = fitCCBehavControlModel(adata,1,strs{si},[],[],lapse);
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
if strfind(modes,'thresholds')
    %% dispInfo
    load(fullfile(datafolder,sprintf('%s',subj),files(end).name));
    ccDispInfocontrol(stimulus,subj);
end
%%
if strfind(modes,'right')
    h1 = cc_rightchoicecontrol(adata, fit, 0);

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