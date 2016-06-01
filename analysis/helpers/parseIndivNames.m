
function [conVal, cohVal, timeVal] = parseIndivNames(conNames,cohNames,timNames)
%% Note that parse names ignores the prefix entirely, so it's hemisphere independent
conVal = []; cohVal = []; timeVal = [];

conStr = 'contrast=';
cohStr = 'coherence=';
timStr = 'timing=';

for i = 1:length(conNames)
    name = conNames{i};
    conVal = [conVal str2num(name(strfind(name,conStr)+length(conStr):end))];
end
for i = 1:length(cohNames)
    name = cohNames{i};
    cohVal = [cohVal str2num(name(strfind(name,cohStr)+length(cohStr):end))];
end
for i = 1:length(timNames)
    name = timNames{i};
    timeVal = [timeVal str2num(name(strfind(name,timStr)+length(timStr):end))];
end