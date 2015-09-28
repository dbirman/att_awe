function fdata = cc_addNeuralFull(allData,name,sid)

data = allData.neural.fit.(name);    

% we are going to convert the data into long form to push it into fullData

% amplitude | task | contrast | coherence | roi
fdata = [];
for ri = 1:length(allData.neural.shortROIs)
    roi = allData.neural.shortROIs{ri};
    [conVals, cohVals, cuedTask] = parseNames(allData.neural.SCM.(name).(roi).stimNames); 
    fdata(:,1) = data.(roi).mat{1}(:);
    fdata(:,2) = cuedTask';
    fdata(:,3) = conVals';
    fdata(:,4) = cohVals';
    fdata(:,5) = repmat(ri,size(cuedTask'));
end





function [conVal, cohVal, cuedTask] = parseNames(stimNames)
%% Note that parse names ignores the prefix entirely, so it's hemisphere independent
conVal = []; cohVal = []; cuedTask = [];

for i = 1:length(stimNames)
    name = stimNames{i};
    ands = strfind(name,' and');
    
    if strfind(name,'Con=')
        conVal(end+1) = str2num(name(strfind(name,'Con=')+4:ands(1))); % note contrast first
    end
    if strfind(name,'Coh=')
        cohVal(end+1) = str2num(name(strfind(name,'Coh=')+4:ands(2)-2));
    end
    % t can be 1, 2, (coherence, contrast) or 3, 4 (catch coherence, catch
    % contrast)
    cuedTask = [cuedTask str2num(name(strfind(name,'ask=')+4:end))];
end

