%%


% subjects = {'s021','s025','s300','s302','s304','s305','s309','s310'};
subjects = {'s0300'};

folders = {{'20150509','20150511','20150513','20150618','20150729'}};
r2_cutoffs = {[0.1,.08,.08,.07,0.08]};
concatScans = {2};
ROIs = {'l_v1','l_v2v','l_v2d','l_v3v','l_v3d','l_hmt','r_v1','r_v2v','r_v2d','r_v3v','r_v3d','r_hmt'};
shortROIs = {'v1','v2','v3','hmt'};

%% Add fields
allData = loadAllData();
for si = 1:length(subjects)
    sid = subjects{si};
    if ~isfield(allData,sid)
        allData.(sid) = struct;
        allData.(sid).neural = struct;
    end
    allData.(sid).neural.folders = folders{si};
    allData.(sid).neural.r2_cutoffs = r2_cutoffs{si};
    allData.(sid).neural.concatScan = concatScans{si};
    allData.(sid).neural.ROIs = ROIs;
    allData.(sid).neural.shortROIs = shortROIs;
    
end
saveAllData(allData);

%% Modify stimfiles
allData = loadAllData();
for si = 1:length(subjects)
    sid = subjects{si};
    
    folders = allData.(sid).neural.folders;
    
    for fi = 1:length(folders)
        fullFolder = fullfile(sprintf('~/data/cohcon/%s%s',sid,folders{fi}));
        disp(sprintf('(matlab_neural) Checking if folder %s needs corrected stimfiles.',folders{fi}));
        genCorrectedStims(fullFolder);
    end
end
% modify doesn't change allData

%% Generate SCM

allData = loadAllData();
for si = 1:length(subjects)
    sid = subjects{si};
    
    allData.(sid).neural.SCM = cc_genSCM(allData.(sid).neural,sid);
end
saveAllData(allData);

%% Simplify SCM

allData = loadAllData();
for si = 1:length(subjects)
    sid = subjects{si};
    
    allData.(sid).neural.SCM = cc_simplifySCM(allData.(sid).neural.SCM,[.2 .4 .6 .8],[0 .02 .1 .2 .4],1,'main');
end
saveAllData(allData);

%% Load ER_analysis

allData = loadAllData();
for si = 1:length(subjects)
    sid = subjects{si};
    
    allData.(sid).neural = cc_loadER(allData.(sid).neural,sid,'main');
end
saveAllData(allData);

%% Concatenate Folders

allData = loadAllData();
for si = 1:length(subjects)
    sid = subjects{si};
    
    allData.(sid).neural = cc_concatER(allData.(sid).neural,'main');
end
saveAllData(allData);

%% Remove stimVols that don't have enough data

allData = loadAllData();
for si = 1:length(subjects)
    sid = subjects{si};
    
    restore = 0; % flag = 1 will reverse the removal and restore the backup
    allData.(sid).neural.SCM = cc_removeNoDataStimVols(allData.(sid).neural.SCM,'main',15,restore);
end
saveAllData(allData);

%% Run Analysis

allData = loadAllData();
for si = 1:length(subjects)
    sid = subjects{si};
    [allData.(sid).neural, allData.(sid).CRF] = cc_performAnalysis(allData.(sid).neural,'main');
end
saveAllData(allData);

%% Save
