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

%% Get all file names
allFiles = {};

for si = 1:length(subjIDs)
    cSub = subjIDs{si};
    for mfi = 1:length(zeroPerf.(cSub).zero_perf.filesM)
        allFiles{end+1} = fullfile('/proj/att_awe/images/brazil_faces/m/',zeroPerf.(cSub).zero_perf.filesM(mfi).name);
    end
    for ffi = 1:length(zeroPerf.(cSub).zero_perf.filesF)
        allFiles{end+1} = fullfile('/proj/att_awe/images/brazil_faces/f/',zeroPerf.(cSub).zero_perf.filesF(ffi).name);
    end
end

%% Remove?
figure
colormap('gray');

for fi = 1:length(allFiles)
    cfile = allFiles{fi};
    img = imread(cfile);
    imagesc(img,[0 255]);
    title(cfile);
    keyboard
    if isequal(input('Remove? ','s'),'s')
        ??
    end
end