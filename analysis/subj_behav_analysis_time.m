function subj_behav_analysis_time( subj, modes )

disp(sprintf('Started running for %s',subj));
%% Get files
files = dir(fullfile(datafolder,'time',subj));

%% Load data
adata = loadadata_time(subj);

%% load original models
if isfile(fullfile(datafolder,'time',sprintf('%s_data.mat',subj)))
    load(fullfile(datafolder,'time',sprintf('%s_data.mat',subj)));
    fits = data.fits;
    BICs = data.BICs;
end

%% Fit Models
if strfind(modes,'refit')
    %% Fit Contrast/Coherence response models (just to control condition)
%     strs = {'con-exp,coh-exp','con-exp,coh-exp,poisson','con-linear,coh-linear','con-linear,coh-linear,poisson'}; %
    strs = {'con-exp,coh-exp','con-exp,coh-exp,timelin','con-exp,coh-exp,notime'};
    %fits = cell(1,length(strs));
%     BICs = zeros(size(fits));
    minl = inf;
    for si = 1:length(strs)
        fits{si} = fitCCBehavControlModel_time(adata,1,strs{si});
        BICs(si) = fits{si}.BIC;
        if fits{si}.BIC < (minl-5)
            minl = fits{si}.BIC;
            fit = fits{si};
        end
    end
end

%% Save data
if exist(fullfile(datafolder,'time',sprintf('%s_data.mat',subj)))==2
    load(fullfile(datafolder,'time',sprintf('%s_data.mat',subj)));
end
if strfind(modes,'refit')
    data.fit = fit;
    data.fits = fits;
    data.BICs = BICs;
    data.strs = strs;
else
    fit = data.fit;
end
save(fullfile(datafolder,'time',sprintf('%s_data.mat',subj)),'data');
% 

%%
disp(sprintf('Finished running for %s',subj));
