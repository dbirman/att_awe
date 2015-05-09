function cohCon_fits2csv(nfits,mfits,cfits)

global analysis
%% This function writes a single file mainThresholds.csv
% This file includes the estimated thresholds for each run of each of the
% task types. The goal is to put these in to LONG form BEFORE exporting to
% CSV, so that the code in R is simplified.

saveLoc = fullfile(analysis.anFolder,'csv');

if ~isdir(saveLoc)
    mkdir(saveLoc);
end

ploFile = fullfile(saveLoc,'modelFits.csv');

%% Plotting(CUES 1/4,PEDESTALS 1:3,DUAL 0/1,NUM noise/con)
% We're going to loop through everything and add it to a matrix that tracks
% the variables in long format instead of in matrix format:

header = {'task','catch','T','beta','lambda','gamma'};
data = [];

for t = 1:2
    for i = 1:size(nfits{t,1},1)
        dat = nfits{t,1}(i,:);
        data(end+1,:) = [t -1 dat(1) dat(2) dat(3) dat(4)];
    end
    for i = 1:size(mfits{t,1},1)
        dat = mfits{t,1}(i,:);
        data(end+1,:) = [t 0 dat(1) dat(2) dat(3) dat(4)];
    end
    for i = 1:size(cfits{t,1},1)
        dat = cfits{t,1}(i,:);
        data(end+1,:) = [t 1 dat(1) dat(2) dat(3) dat(4)];
    end
end

%% Write file

csvwriteh(ploFile,data,header);