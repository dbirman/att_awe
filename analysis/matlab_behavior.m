%% Setup


% subjects = {'s021','s025','s300','s302','s304','s305','s309','s310'};
subjects = {'s304','s305','s310','s315'};

%% Copy Files

for si = 1:length(subjects)
    sid = subjects{si};
    
    cc_copyFiles(sid);
end

%% Load and Export

allData = loadAllData(sid);

maf = {};
caf = {};
nocaf = {};
for si = 1:length(subjects)
    sid = subjects{si};
    
    mglSetSID(sid);
    
    allData = exportCohCon(allData);
    saveAllData(sid,allData);
end


%% Eye Plots

for si = 1:length(subjects)
    sid = subjects{si};
    
    cc_eyePlots(sid);
end

%% R-Choice Plots

for si = 1:length(subjects)
    sid = subjects{si};
    
    cc_rightChoicePlots(sid, false);   
end

%% State-Trace Plots

for si = 1:length(subjects)
    sid = subjects{si};
        
    [f_stp.(sid).f, f_stp.(sid).data] = cc_stateTracePlot(sid);
end

%% Simulations

simData = cc_runBehavSimulations();

%% Generate fullData
% fullData is the entire dataset in long form, it loads all of the subjects
% individually and copies their entire dataset (every trial) into a single
% matrix. The formatting is as follows:
% Run # | Trial # | lCon | rCon | lCoh | rCoh | dCon | dCoh | side | response | task 

fullData = loadFullData();

if length(fields(fullData))==0 || isempty(fullData.data) || ~length(subjects)==length(fullData.subjects)
    fullData.subjects = subjects;
    fData = {};
    for si = 1:length(subjects)
        sid = subjects{si};

        [fData{si}, header] = cc_exportFull(sid);
        fData{si}(:,end+1) = mrStr2num(strrep(sid,'s',''));
    end
    header{end+1} = 'sid';

    data = cc_catLong(fData);

    % add 'side'
    header{end+1} = 'side';

    side = [];
    for i = 1:size(data)
        if data(i,3)==1 % coherence
            side(i) = data(i,10);
        else
            side(i) = data(i,9);
        end
    end
    data(:,end+1) = side;

    fullData.data = data;
    fullData.header = header;
end

saveFullData(fullData);

%% lightning talk plot

colors = brewermap(15,'PuOr');

mafs = []; cafs = [];
for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    
    figure
    subplot(1,2,1), hold on
    % attending coherence
    bar([1 2],[allData.behav.maf{1}.threshold allData.behav.caf{1}.threshold]);
    axis([0.5 2.5 0 .4]);
    set(gca,'yscale','log');
    subplot(1,2,2), hold on
    % attending contrast
    bar([1 2],[allData.behav.maf{2}.threshold allData.behav.caf{2}.threshold]);
    axis([0.5 2.5 0 0.5]);
    set(gca,'yscale','log');
    dir = fullfile('~/proj/att_awe/analysis/figures',sid);
    if ~isdir(dir), mkdir(dir); end
    fname = fullfile(dir,'threshold_bar.pdf');
    print(fname,'-dpdf');
end


%% junk

% % % % % 
% % % % % %% Generate CDF for selectionModel
% % % % % %        CDF is a similar struncture containing contrast discrimination data. So, continuing
% % % % % %        from the exmaple above, if you have a focal cue condition you would have fields for
% % % % % %        both data and fit that contain contrast discrimination data:
% % % % % %        - CDF.focal.data.c     = [0.01 0.5 0.7];  % Pedestal contrast 
% % % % % %        - CDF.focal.data.t     = [0.04 0.08 0.12];% Threshold, delta-c
% % % % % %        - CDF.focal.data.tste  = [0.01 0.02 0.01];% ste of threshold
% % % % % %        - CDF.focal.data.d     = {[0.0 0.01 0.02], ...
% % % % % %                                  [0.01 0.5 0.7],  ...
% % % % % %                                  [0.01 0.5 0.7]};
% % % % % %        - CDF.focal.data.nCued = 1;
% % % % % % 
% % % % % %        Where c are the pedestal contrasts, t are the measured thresholds and d is the
% % % % % %        distractor contrasts in a cell array. The cued field contains how many targets
% % % % % %        were cued. The tste field contains the standard error of the CDF data (for display)
% % % % % for task = 1:2
% % % % %     taskStr = {'coherence','contrast'};
% % % % %     peds = [.1 .6];
% % % % %     
% % % % %     % we copy for each roi
% % % % %     CDF.(taskStr{task}).nocatch.pedestal = peds(task);
% % % % %     CDF.(taskStr{task}).nocatch.t = maf{task}.fit.fitparams(1);
% % % % %     CDF.(taskStr{task}).nocatch.delta = maf{task}.fit.signal;
% % % % %     CDF.(taskStr{task}).nocatch.p = maf{task}.fit.pcorrect;
% % % % %     CDF.(taskStr{task}).nocatch.pste = maf{task}.fit.pcorrectste;
% % % % %     
% % % % % end