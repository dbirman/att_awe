%% Analysis v2.0
% Okay, I've gone through a number of ideas and have settled on two basic
% analyses to look at the data. The idea is to build both an encoding
% model, which we will use as a general analysis of what brain areas are
% involved in the task in each attention condition, as well as a decoding
% model which we will use a specific test of the ability to the task.

%% Encoding Model
% For a timeseries, in a given attention condition, build the following model:
% Amplitude ~ Contrast + Coherence + Resp
%
% By building this model at three levels:
% Individual Voxels
% Searchlight (across flatmap)
% ROI
%
% We will be able to analyze in a non-ROI based manner how the brain is
% encoding the different parts of the task in each of the attentional
% conditions.
%
% Each of these three 

%% Decoding Model

% For each attention condition build the following models:
% Contrast ~ Voxel Amplitudes (across ROI)
% Coherence ~ Voxel Amplitudes (across ROI)
%
% We will use these models to predict behavior in correct trials in the
% following way:
% (1) Collect instances x voxels for several sizes:
%   a - by ROI
%   b - across hemispheres
% (2) Leave-one-out cross validation