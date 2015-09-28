%%


asubjects = {'s0300','s0304','s0315','s0305'};

afolders = {{'20150509','20150511','20150513','20150618','20150729'},...
    {'20150729','20150801','20150804a','20150804b','20150810'},...
    {'20150906','20150907','20150908','20150909'},...
    {'20150815','20150818','20150822','20150830'}};
ar2_cutoffs = {[0.12,.1,.1,.09,0.1],...
    [0.15,0.08,0.1,0.1,0.09],...
    [0.1,0.15,0.1,0.15,0.1],...
    [0.07,0.05,0.07,0.08]};
aconcatScans = {2,1,1,1};

ROIs = {'lV1','lV2v','lV2d','lV3v','lV3d','lMT','rV1','rV2v','rV2d','rV3v','rV3d','rMT','lV3a','lV3b','rV3a','rV3b'};
shortROIs = {'V3a','V3b','V1','V2','V3','MT'};

reset = [0 0 0 1];

run = [4];

%% Use run to alter fields

subjects = asubjects(run);
folders = afolders(run);
r2_cutoffs = ar2_cutoffs(run);
concatScans = aconcatScans(run);
reset = reset(run);

%% Add fields
for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    if reset(si)
        allData = struct;
    end
    if ~isfield(allData,'neural')
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

for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    allData.neural = cc_neuralFigures(allData.neural,'highr2',sid,false,true);
    saveAllData(sid,allData);
end

%% Full Data
disp('pausing');
keyboard
% here we want to save the neural data points for the full model analysis.
% What we do need? We'll keep the deconvolved response to each of the
% binned conditions, so we end up with like a 4x4 matrix of deconvolved
% responses

header = {'amplitude','task','contrast','coherence','roi','sid'};
slong = {};
for si = 1:length(asubjects)
    sid = asubjects{si};
    allData = loadAllData(sid);
    slong{si} = cc_addNeuralFull(allData,'highr2',sid);
end
fullN = cc_catLong(slong); fullN.header = header;

%% Cat full
fullData = loadFullData;
fullData.ndata = fullN.data;
fullData.nheader = fullN.header;
fullData.nsubjects = asubjects;
saveFullData(fullData);