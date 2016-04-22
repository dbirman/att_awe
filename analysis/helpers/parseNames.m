
function [conVal, cohVal] = parseNames(stimNames,constr,cohstr,sep)
%% Note that parse names ignores the prefix entirely, so it's hemisphere independent
conVal = []; cohVal = []; cuedTask = [];

for i = 1:length(stimNames)
    name = stimNames{i};
    ands = strfind(name,sep);
    
    if ~isempty(strfind(name,constr))
        conVal(end+1) = str2num(name(strfind(name,constr)+length(constr):ands(1))); % note contrast first
    end
    if ~isempty(strfind(name,cohstr))
        cohVal(end+1) = str2num(name(strfind(name,cohstr)+length(cohstr):end));
    end
end
