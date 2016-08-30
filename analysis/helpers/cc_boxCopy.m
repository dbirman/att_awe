%% Definitions
% subject lists
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];
% folder locations
nData = '~/data/cohcon_localizer';
nSuffix = '%s*_data.mat';
bData = '~/data/cohcon';
bSuffix = '%s/*_stim*.mat';
boxsync = '/Users/dan/proj/Box Sync/COHCON_DATA';

% This script copies files into the Box Sync folder that is shared between
% your MAC and PC computer.

%% Copy Behavior Files
% this copies every single stim file
for ni = 1:length(aSIDs)
    sid = aSIDs(ni);
    subj = sprintf('s%03.0f',sid);
    disp(sprintf('### SUBJ: %s ###',subj));
    cc_copyFiles(subj);
    files = dir(fullfile(bData,sprintf(bSuffix,subj)));
    copyDir = fullfile(boxsync,subj);
    if ~isdir(copyDir), mkdir(copyDir); end
    for fi = 1:length(files)
        fname = fullfile(bData,subj,files(fi).name);
        if ~isfile(fullfile(copyDir,files(fi).name))
            copyfile(fname,fullfile(copyDir,files(fi).name));
        end
    end
end

%% Run savedata_localizers
version = 2;
for ni = 1:length(nSIDs)
    % get folders
    subj = nSIDs(ni);
    files = dir(sprintf('~/data/cohcon_localizer/s%04.0f*',subj));
    for fi = 1:length(files)
        disp(sprintf('Starting %s',files(fi).name));
        if isempty(strfind(files(fi).name,'mat')) && isempty(strfind(files(fi).name,'_'))
            skip = false;
            if isfile(fullfile(nData,sprintf('%s_data.mat',files(fi).name)))
                load(fullfile(nData,sprintf('%s_data.mat',files(fi).name)))
                if isfield(data,'version') && data.version==version;
                    skip = true;
                end
            end
            if ~skip
                savedata_localizer(files(fi).name,version);
            else
                disp('Skipping');
            end
        end
    end
end

%% Copy Neural Files
for ni = 1:length(nSIDs)
    sid = nSIDs(ni);
    subj = sprintf('s%04.0f',sid);
    disp(sprintf('### SUBJ: %s ###',subj));
    files = dir(fullfile(nData,sprintf(nSuffix,subj)));
    for fi = 1:length(files)
        fname = fullfile(nData,files(fi).name);
        copyfile(fname,fullfile(boxsync,files(fi).name));
    end
end