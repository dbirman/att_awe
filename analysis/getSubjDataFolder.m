function [datFolder, anFolder] = getSubjDataFolder( )

if isempty(mglGetSID)
    error('Please set subject ID first.');
end

datFolder = sprintf('~/data/cue_noisecon/%s',mglGetSID);
anFolder = sprintf('~/proj/att_awe/analysis/%s',mglGetSID);

if ~isdir(anFolder)
    mkdir(anFolder);
end
if ~isdir(fullfile(anFolder,'figures'))
    mkdir(fullfile(anFolder,'figures'));
end
