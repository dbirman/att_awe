%% Setup


% subjects = {'s021','s025','s300','s302','s304','s305','s309','s310'};
subjects = {'s304'};

%% Copy Files

for si = 1:length(subjects)
    sid = subjects{si};
    
    cc_copyFiles(sid);
end

%% Load and Export

allData = loadAllData();

maf = {};
caf = {};
nocaf = {};
for si = 1:length(subjects)
    sid = subjects{si};
    
    mglSetSID(sid);
    
    allData = exportCohCon(allData);
end

saveAllData(allData);

%% R-Choice Plots

cc_rightChoicePlots(allData,sid);

%% State-Trace Plots

for si = 1:length(subjects)
    sid = subjects{si};
        
    [f_stp.(sid).f, f_stp.(sid).data] = cc_stateTracePlot(allData,sid);
end

%% Simulations

simData = cc_runBehavSimulations();
%%

%% Generate CDF for selectionModel
%        CDF is a similar struncture containing contrast discrimination data. So, continuing
%        from the exmaple above, if you have a focal cue condition you would have fields for
%        both data and fit that contain contrast discrimination data:
%        - CDF.focal.data.c     = [0.01 0.5 0.7];  % Pedestal contrast 
%        - CDF.focal.data.t     = [0.04 0.08 0.12];% Threshold, delta-c
%        - CDF.focal.data.tste  = [0.01 0.02 0.01];% ste of threshold
%        - CDF.focal.data.d     = {[0.0 0.01 0.02], ...
%                                  [0.01 0.5 0.7],  ...
%                                  [0.01 0.5 0.7]};
%        - CDF.focal.data.nCued = 1;
% 
%        Where c are the pedestal contrasts, t are the measured thresholds and d is the
%        distractor contrasts in a cell array. The cued field contains how many targets
%        were cued. The tste field contains the standard error of the CDF data (for display)
for task = 1:2
    taskStr = {'coherence','contrast'};
    peds = [.1 .6];
    
    % we copy for each roi
    CDF.(taskStr{task}).nocatch.pedestal = peds(task);
    CDF.(taskStr{task}).nocatch.t = maf{task}.fit.fitparams(1);
    CDF.(taskStr{task}).nocatch.delta = maf{task}.fit.signal;
    CDF.(taskStr{task}).nocatch.p = maf{task}.fit.pcorrect;
    CDF.(taskStr{task}).nocatch.pste = maf{task}.fit.pcorrectste;
    
end