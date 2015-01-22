function plo2csv(plotting)

%% This function writes a single file mainThresholds.csv
% This file includes the estimated thresholds for each run of each of the
% task types. The goal is to put these in to LONG form BEFORE exporting to
% CSV, so that the code in R is simplified.

saveLoc = fullfile('~/proj/att_awe/analysis/',mglGetSID,'csv');

if ~isdir(saveLoc)
    mkdir(saveLoc);
end

ploFile = fullfile(saveLoc,'mainThresholds.csv');

%% Plotting(CUES 1/4,PEDESTALS 1:3,DUAL 0/1,NUM noise/con)
% We're going to loop through everything and add it to a matrix that tracks
% the variables in long format instead of in matrix format:

header = {'estimate','task','cues','pedestal','dual','threshold'};
data = [];

for c = 1:2
    for p = 1:3
        for d = 1:2
            for t = 1:2
                threshs = plotting{c,p,d,t};
                for i = 1:length(threshs)
                    data(end+1,:) = [i t c p d threshs(i)];
                end
            end
        end
    end
end

%% Write file

csvwriteh(ploFile,data,header);