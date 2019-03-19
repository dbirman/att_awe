
%% Helper script for copying data:

nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];

% copy all the files
for ai = 1:length(aSIDs)
    copyfile(fullfile(datafolder,sprintf('s%03.0f_adata.mat',aSIDs(ai))),fullfile('~/proj/ej_class/data/',sprintf('%i_data.mat',ai)));
end