function fullData = loadFullData( )

fFile = fullfile('~/proj/att_awe/analysis/data/fullData.mat');

if isfile(fFile)
    disp('(fullData) Loading...');
    load(fFile);
else
    disp('(fullData) Generating New...');
    fullData = struct;
end