function [fullData, header] = cc_exportFull(sid)

mglSetSID(sid);

doEye = false;
scan = false;

global analysis

[analysis.datFolder, analysis.anFolder] = getSubjDataFolder('cohcon',scan);

year = date;
year = year(end-1:end);

files = dir(sprintf('%s/%s*mat',analysis.datFolder,year));

%% Loop over files and load them into expHolder

expHolder = loadExp(files);


%% Loop over exp files and export to a single .mat file
fData = {};
for ei = 1:length(expHolder)
    [fData{ei} header] = cc_2long(expHolder{ei});
end
fullData = cc_catLong(fData);