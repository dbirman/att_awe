%% Copy All Data + Analysis Scripts

% This script copies all the data and analysis scripts in the
% ~/data/cohcon_localizer folder into ~/proj/att_awe so that they live on
% github and can be analyzed separately.

afiles = dir('~/data/cohcon_localizer/*.m');
mfiles = dir('~/data/cohcon_localizer/*.mat');

%%
if ~isdir('~/proj/att_awe/model'), mkdir('~/proj/att_awe/model'); end

%%
for ai = 1:length(afiles)
    copyfile(fullfile('~/data/cohcon_localizer',afiles(ai).name),fullfile('~/proj/att_awe/model',afiles(ai).name));
end

for mi = 1:length(mfiles)
    copyfile(fullfile('~/data/cohcon_localizer',mfiles(mi).name),fullfile('~/proj/att_awe/model',mfiles(mi).name));
end