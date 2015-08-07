function allData = loadAllData( )
%LOADALLDATA Summary of this function goes here
%   Detailed explanation goes here

adFile = '~/proj/att_awe/analysis/allData.mat';

if isfile(adFile)
    disp('(allData) Loading...');
    load(adFile);
else
    disp('(allData) Generating New...');
    allData = struct;
end

end

