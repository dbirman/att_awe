function allData = loadAllData(sid)
%LOADALLDATA Summary of this function goes here
%   Detailed explanation goes here

adFile = fullfile('~/proj/att_awe/analysis/data/',sprintf('%s_allData.mat',sid));

if isfile(adFile)
    disp('(allData) Loading...');
    load(adFile);
else
    disp('(allData) Generating New...');
    allData = struct;
end

end
