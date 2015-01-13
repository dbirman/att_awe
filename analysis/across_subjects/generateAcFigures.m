% generateAcFigures
%   Generates across subject figures for analysis

%% Load zero performance images data

subjIDs = {'s025','s300'};

zeroPerf = struct;
for si = 1:length(subjIDs)
    cSub = subjIDs{si};
    mglSetSID(cSub);
    [datFolder, anFolder] = getSubjDataFolder;
    zeroPerf.(cSub) = load(fullfile(anFolder,'zero_perf.mat'));
end

clear datFolder anFolder

%% Check for overlap

for si = 1:length(subjIDs)
    cSub = subjIDs{si};