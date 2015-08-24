%%


% subjects = {'s021','s025','s300','s302','s304','s305','s309','s310'};
asubjects = {'s0300','s0304'};

folders = {{'20150509','20150511','20150513','20150618','20150729'},...
    {'20150729','20150801','20150804a','20150804b','20150810'}};
r2_cutoffs = {[0.12,.1,.1,.09,0.1],...
    [0.15,0.08,0.1,0.1,0.09]};
concatScans = {2,1};

%%%%% TODO %%%%% EDIT ROIs to reflect naming conventions
ROIs = {'l_v1','l_v2v','l_v2d','l_v3v','l_v3d','l_hmt','r_v1','r_v2v','r_v2d','r_v3v','r_v3d','r_hmt'};
shortROIs = {'v1','v2','v3','hmt'};

run = [2];

%% Use run to alter fields

subjects = asubjects(run);
folders = folders(run);
r2_cutoffs = r2_cutoffs(run);
concatScans = concatScans(run);

%% Add fields
for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    if ~isfield(allData,sid)
        allData = struct;
        allData.neural = struct;
    end
    allData.neural.folders = folders{si};
    allData.neural.r2_cutoffs = r2_cutoffs{si};
    allData.neural.concatScan = concatScans{si};
    allData.neural.ROIs = ROIs;
    allData.neural.shortROIs = shortROIs;
    saveAllData(sid,allData);
end

%% Modify stimfiles
for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    
    folders = allData.neural.folders;
    
    for fi = 1:length(folders)
        fullFolder = fullfile(sprintf('~/data/cohcon/%s%s',sid,folders{fi}));
        disp(sprintf('(matlab_neural) Checking if folder %s needs corrected stimfiles.',folders{fi}));
        genCorrectedStims(fullFolder);
    end
end
% modify doesn't change allData

%% Generate SCM

for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    
    disp('(Rewrite here: copy existing SCM if it exists)');
    allData.neural.SCM = cc_genSCM(allData.neural,sid);
    saveAllData(sid,allData);
end

%% Simplify SCM

for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    
    allData.neural.SCM = cc_simplifySCM(allData.neural.SCM,[.2 .4 .6 .8 1],[0 .02 .1 .2 .4 .6],1,'highr2');
    saveAllData(sid,allData);
end
%% Load ER_analysis

for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    
    allData.neural = cc_loadER(allData.neural,sid,'highr2');
    saveAllData(sid,allData);
end

%% Concatenate Folders

for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    
    allData.neural = cc_concatER(allData.neural,'highr2');
    saveAllData(sid,allData);
end
%% Remove stimVols that don't have enough data

for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    
    restore = 0; % flag = 1 will reverse the removal and restore the backup
    allData.neural.SCM = cc_removeNoDataStimVols(allData.neural.SCM,'highr2',15,restore);
    saveAllData(sid,allData);
end
%% Run Analysis

for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    [allData.neural] = cc_performAnalysis(allData.neural,'highr2');
    saveAllData(sid,allData);
end

%% Figures

for si = 1:length(asubjects)
    sid = asubjects{si};
    allData = loadAllData(sid);
    [allData.CRF] = cc_neuralFigures(allData.neural,'highr2');
    saveAllData(sid,allData);
end