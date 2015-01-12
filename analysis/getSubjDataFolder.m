function datFolder = getSubjDataFolder( )

if isempty(mglGetSID)
    error('Please set subject ID first.');
end

datFolder = sprintf('~/data/cue_noisecon/%s',mglGetSID);
