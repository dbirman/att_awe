function fade_cat2csv(peripheral)
global analysis

%% This function writes a single file perThresholds.csv
% This file includes the estimated thresholds for each run of each of the
% task types. The goal is to put these in to LONG form BEFORE exporting to
% CSV, so that the code in R is simplified.

saveLoc = fullfile(analysis.anFolder,'csv');

if ~isdir(saveLoc)
    mkdir(saveLoc);
end

perFile = fullfile(saveLoc,'catThresholds.csv');

%% Peripheral(TASK single, dual:noise, dual:contrast)
% We're going to loop through everything and add it to a matrix that tracks
% the variables in long format instead of in matrix format:

header = {'estimate','task','threshold'};
data = [];

for t = 1:3
    cur = peripheral{t};
    for i = 1:length(cur)
        data(end+1,:) = [i t cur(i)];
    end
end

%% Write file

csvwriteh(perFile,data,header);