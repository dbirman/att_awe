
function [conVal, cohVal, timeVal] = parseNames(stimNames,constr,cohstr,timestr,sep)
%% Note that parse names ignores the prefix entirely, so it's hemisphere independent
conVal = []; cohVal = []; timeVal = [];

for i = 1:length(stimNames)
    name = stimNames{i};
    if ~isempty(sep)
        ands = strfind(name,sep);
    else
        ands = length(name);
    end
    
    andused = 0;
    
    % we always do contrast first when possible
    if ~isempty(constr)
        if ~isempty(strfind(name,constr))
            conVal(end+1) = str2num(name(strfind(name,constr)+length(constr):ands(1))); % note contrast first
            andused = 1;
        end
    end
    % coherence comes next
    if ~isempty(cohstr)
        if ~isempty(strfind(name,cohstr))
            if andused
                cohVal(end+1) = str2num(name(strfind(name,cohstr)+length(cohstr):end));
            else
                cohVal(end+1) = str2num(name(strfind(name,cohstr)+length(cohstr):ands(1)));
                andused = 1;
            end
        end
    end
    % timing is always last
    if ~isempty(timestr)
        if ~isempty(strfind(name,timestr))
            timeVal(end+1) = str2num(name(strfind(name,timestr)+length(timestr):end));
        end
    end
end
