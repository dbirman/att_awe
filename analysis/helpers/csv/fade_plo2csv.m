function fade_plo2csv(plotting)

global analysis
%% This function writes a single file mainThresholds.csv
% This file includes the estimated thresholds for each run of each of the
% task types. The goal is to put these in to LONG form BEFORE exporting to
% CSV, so that the code in R is simplified.

saveLoc = fullfile(analysis.anFolder,'csv');

if ~isdir(saveLoc)
    mkdir(saveLoc);
end

ploFile = fullfile(saveLoc,'mainThresholds.csv');

%% Plotting(CUES 1/4,PEDESTALS 1:3,DUAL 0/1,NUM noise/con)
% We're going to loop through everything and add it to a matrix that tracks
% the variables in long format instead of in matrix format:

header = {'task','cues','pedestal','threshold'};
data = [];

t = 1; % always contrast
for c = 1:2
    for p = 1:3
        threshs = plotting{c,p};
        for i = 1:length(threshs)
            data(end+1,:) = [t c p threshs(i)];
        end
    end
end

%% Write file

csvwriteh(ploFile,data,header);