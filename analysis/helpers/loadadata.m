function adata = loadadata(subj)
% format:
%     1       2         3        4      5      6     7      8       9
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch -
%      10      11        12
%   pedcon - pedcoh - correct

if strfind(getenv('OS'),'Windows')
    cfolder = sprintf('C:/Users/Dan/proj/COHCON_DATA/%s',subj);
else
    cfolder = sprintf('~/proj/data/cohcon/%s',subj);
end

files = dir(fullfile(cfolder,'*.mat'));

%% 
adata = [];
for fi = 1:length(files)
    load(fullfile(cfolder,files(fi).name));
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
disp(sprintf('%s trials %i',subj,size(adata,1)));

%% Remove NaN
adata = adata(~any(isnan(adata),2),:);