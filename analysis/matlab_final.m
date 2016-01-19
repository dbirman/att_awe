
%% Load allData
adata = loadAllData(sid);
neural = adata.neural;

%% Analysis v2.0
% Okay, I've gone through a number of ideas and have settled on two basic
% analyses to look at the data. The idea is to build both an encoding
% model, which we will use as a general analysis of what brain areas are
% involved in the task in each attention condition, as well as a decoding
% model which we will use a specific test of the ability to the task.


%% Decoding Model
% do this analysis first

% For each attention condition build the following models:
% Contrast ~ Voxel Amplitudes (across ROI)
% Coherence ~ Voxel Amplitudes (across ROI)
% Resp ~ Voxel Amplitudes (across brain?? not sure how to do this)
%
% We will use these models to predict behavior in correct trials in the
% following way:
% (1) Collect instances x voxels for several sizes:
%   a - by ROI
%   b - across hemispheres
% (2) Leave-one-out cross validation
%   predict each trials outcome from the difference in estimated contrast
%   from left/right hemispheres
% (3) Create a Sensitivity Map
%   fit full model (not leave-one-out) and then use the beta values for 
%   each hemispheric ROI, to show that different areas are sensitive in the
%   different attention conditions 
%
% What do we need:
% ROIs loaded for each hemisphere
% stimvols, amp/con/coh
%
% What do we do:
% train CV lasso and pick best value, use to predict L1O, continue
% iterating across examples
%
% Output:
% (1)
%   performance across instances for each attentional condition (Which should
%   map onto the behavioral data for that participant)
% (2)
%   a map of sensitivity values across cortex, showing that V1/MT are
%   sensitive in different tasks
%
%
%

%% Decode: Preload Data for Subject + Stimvols (convert to .mat)
ds = cc_load(neural,sid);

%% Decode: Organize Data into ROI data files
insts = cc_inst(neural,sid,ds);

%% Decode: Build the lasso dataset
datasets = cc_buildset(neural,sid,insts,{'V1','V2','V3','MT','V4'});

%% combine datasets
fldata = []; frdata = [];
for fi = 1:length(datasets)
    load(datasets{fi});
    %ldata/rdata
    fldata = [fldata;ldata];
    frdata = [frdata;rdata];
end
ldata = fldata; rdata = frdata;
save('~/data/cohcon/s0300_dataset_concat.mat','ldata','rdata');
clear fldata frdata rdata ldata
dsconcat = {'~/data/cohcon/s0300_dataset_concat.mat'};

%% Cleanup datasets (remove amplitudes >/< 10/-10)
clean_ds = cc_clean(dsconcat);

%% Decode: Build lasso models
output = cc_lasso(neural,sid,clean_ds);
% procedure is: put test group aside, CV lasso parameters, estimate test
% performance, CV 10x (inside parfor)
% output is: performance, for each folder, for each condition (task=1/2 on
% con/coh estimation)

%% Decode: Produce output maps and send to joint file

%% Encoding Model
% do this analysis second

% For a timeseries, in a given attention condition, build the following model:
% Amplitude ~ Contrast + Coherence + Resp
%
% By building this model at three levels:
% Individual Voxels
% Individual Voxels after a Gaussian Kernel Smoothing (in surface space)
% ROI
%
% We will be able to analyze in a non-ROI based manner how the brain is
% encoding the different parts of the task in each of the attentional
% conditions.
%
% Each of these three analyses create a sensitivity map for contrast,
% coherence, and responses
