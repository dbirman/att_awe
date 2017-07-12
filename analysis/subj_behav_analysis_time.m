function subj_behav_analysis_time( subj, modes )

disp(sprintf('Started running for %s',subj));
%% Get files
files = dir(fullfile(datafolder,'time',subj));

%% Load data
adata = loadadata_time(subj);


%% Fit Models
if strfind(modes,'refit')

    fits{si} = fitCCBehavControlModel_time(adata,0,'freeze',squeeze(respcon(tcorrespond();
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
