function [adata, n] = loadadata(subj)
% format:
%     1       2         3        4      5      6     7      8       9
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch -
%      10      11        12
%   pedcon - pedcoh - correct

files = dir(fullfile(datafolder,subj,'*.mat'));

%% 
adata = [];
for fi = 1:length(files)
    load(fullfile(datafolder,subj,files(fi).name));
    e = getTaskParameters(myscreen,task);
    if length(e{1})==2
        e = e{1}(2);
        flip = [-1 1];
        adata = [adata ; [repmat(stimulus.runs.curTask*flip(stimulus.nocatch+1),length(e.response),1), ...
            repmat(stimulus.baseCon,length(e.response),1), repmat(stimulus.baseCoh,length(e.response),1), ...
            e.randVars.lCon', e.randVars.rCon',...
            e.randVars.lCoh', e.randVars.rCoh',...
            (e.response-1)', (e.parameter.catch)',...
            e.randVars.contrast', e.randVars.coherence',e.randVars.correct']];
    end
end
n = sum(adata(:,9)==-1);
disp(n);
disp(sprintf('%s trials %i',subj,size(adata,1)));

%% Remove NaN
adata = adata(~any(isnan(adata),2),:);