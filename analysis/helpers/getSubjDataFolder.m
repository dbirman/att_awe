function [datFolder, anFolder] = getSubjDataFolder( experiment )

if isempty(mglGetSID)
    error('Please set subject ID first.');
end

switch experiment
    case 'noisecon'
        datFolder = fullfile(sprintf('~/data/cue_noisecon/%s',mglGetSID));
    case 'fade'
        datFolder = fullfile(sprintf('~/data/cue_fade/%s',mglGetSID));
end
        
anFolder = fullfile(datFolder,'analysis');

if ~isdir(anFolder)
    mkdir(anFolder);
end
if ~isdir(fullfile(anFolder,'figures'))
    mkdir(fullfile(anFolder,'figures'));
end
