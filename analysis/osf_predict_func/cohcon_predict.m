function out = cohcon_predict(areas,dcon,dcoh,duration,varargin)
%COHCON_PREDICT Contrast, coherence, and duration predictions for retinotopic
%visual cortex
%
%   Usage: cohcon_predict(areas, delta contrast, delta coherence, duration)
%   By: Daniel Birman, Justin Gardner
%   Date: 4/4/18
%   e.g.: cohcon_predict('V1',0.25,0.5,1);
%   e.g.: cohcon_predict({'V1','V2'},0.25,0.5,1);
%   e.g.: cohcon_predict({'V3A','MT'},0.25,0.5,1,'baseline=[0.5,0]','subjects=1');
%
%   Purpose: Predict the cortical response to a change in contrast and
%       coherence. Responses are in % signal change (fMRI). See paper for
%       additional details.
%
%   Input arguments:
%       Areas: V1, V2, V3, V4, V3A, V3B, V7, MT
%       Delta contrast: 0->1 (can be an array)
%       Delta coherence: 0->1 (can be an array)
%       Duration: 0->10 s (data based on 0.25 to 4 s)
%       subjects (default=0): Set to 1 to get individual subject data, set
%           to an array to get averages
%
%   Output:
%       .contrastResponse: the % signal evoked by a change from baseline to
%           the requested delta contrasts. Organzed as a contrast * area
%           matrix.
%       .coherenceResponse
%       .contrast: copy of inputs
%       .coherence: copy of inputs
%       .visualAreas: copy of inputs
%
%   Please cite:
%       Birman, D., & Gardner, J. L. (2018). A quantitative framework for 
%       motion visibility in human cortex. Journal of neurophysiology.

warning('This function isn''t quite complete -- I still need to add the duration and have it output not only the amplitudes of the response but the convolved hemodynamic response at all. Sorry for the delay -- email dbirman@stanford.edu if you have any questions in the meantime.');
subjects=1:11;
getArgs(varargin,{'baseline=[0.25 0]','subjects=[1:11]'});

if ~isfile('params.mat')
    error('Please download the parameters file params.mat from https://osf.io/s7j9p/');
end

if ~iscell(areas)
    areas = {areas};
end

load('params.mat');

conout = zeros(length(subjects),length(areas),length(dcon));
cohout = zeros(length(subjects),length(areas),length(dcoh));

for ni = subjects
    sparam = params{ni};
    
    for ai = 1:length(areas)
        param = sparam.(areas{ai});
        
        conout(ni,ai,:) = conResp(dcon,param.conalpha,param.consigma,param.conp,param.conq);
        cohout(ni,ai,:) = cohResp(dcoh,param.cohalpha,param.cohkappa);
    end
end

out.contrast = dcon';
out.contrastResponse = squeeze(mean(conout))';
out.coherence = dcoh';
out.coherenceResponse = squeeze(mean(cohout))';
out.visualAreas = areas;

function out = conResp(con,alpha,sigma,p,q)

out = alpha .* ((con.^(p+q)) ./ (con.^q + sigma.^q));

function out = cohResp(coh,alpha,kappa)

out = alpha-(alpha * exp(-kappa*coh));
