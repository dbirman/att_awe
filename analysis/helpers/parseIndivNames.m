
function [conVal, cohVal, timeVal, lConVal, rConVal, lCohVal, rCohVal, taskVal] = parseIndivNames(conNames,cohNames,timNames,lConNames,rConNames,lCohNames,rCohNames,taskNames)
%% Note that parse names ignores the prefix entirely, so it's hemisphere independent
conVal = []; cohVal = []; timeVal = []; lConVal = []; rConVal = []; lCohVal = []; rCohVal = []; taskVal = [];

conStr = 'contrast=';
cohStr = 'coherence=';
timStr = 'timing=';
lconStr = 'lCon=';
rconStr = 'rCon=';
lcohStr = 'lCoh=';
rcohStr = 'rCoh=';
taskStr = 'task=';

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
for i = 1:length(lConNames)
    name = lConNames{i};
    lConVal(end+1) = str2num(name(strfind(name,lconStr)+length(lconStr):end));
end
for i = 1:length(rConNames)
    name = rConNames{i};
    rConVal(end+1) = str2num(name(strfind(name,rconStr)+length(rconStr):end));
end
for i = 1:length(lCohNames)
    name = lCohNames{i};
    lCohVal(end+1) = str2num(name(strfind(name,lcohStr)+length(lcohStr):end));
end
for i = 1:length(rCohNames)
    name = rCohNames{i};
    rCohVal(end+1) = str2num(name(strfind(name,rcohStr)+length(rcohStr):end));
end
for i = 1:length(taskNames)
    name = taskNames{i};
    taskVal(end+1) = str2num(name(strfind(name,taskStr)+length(taskStr):end));
end