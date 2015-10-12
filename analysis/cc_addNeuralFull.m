function fdata = cc_addNeuralFull(allData,name,sid)

data = allData.neural.(name);    

% we are going to convert the data into long form to push it into fullData

% amplitude | se | task | contrast | coherence | roi | sid
fdata = [];
for ri = 1:length(allData.neural.shortROIs)
    ndata = [];
    roi = allData.neural.shortROIs{ri};
    [conVals, cohVals, cuedTask] = parseNames(allData.neural.SCM.(name).(roi).stimNames); 
    ndata(:,1) = data.fits{ri}.amplitude;
    ndata(:,2) = data.fits{ri}.amplitudeSTE;
    ndata(:,3) = cuedTask';
    ndata(:,4) = conVals';
    ndata(:,5) = cohVals';
    ndata(:,6) = repmat(ri,size(cuedTask'));
    ndata(:,7) = repmat(sid,size(cuedTask'));
    fdata = [fdata ; ndata];
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

