function [neural] = cc_performAnalysis(neural,name)

neural.fit = struct;

rois = neural.shortROIs;

if ~isdir(fullfile('~/proj/att_awe/analysis/figures/fitCheck/')), mkdir(fullfile('~/proj/att_awe/analysis/figures/fitCheck/')); end

%% Calculate Non-Lin Fits (for amplitude)

fitter = 'fitType=nonlin';
amper = 'amplitudeType=fit2';
fits = {};
for ri = 1:length(rois)
    roi = rois{ri};       
    
    fit = fitTimecourse(neural.tSeries.(name).(roi).tSeries,neural.SCM.(name).(roi).stimVol,.5,'concatInfo',neural.tSeries.(name).(roi).concatInfo,fitter,amper);
    figure(1)
    fname = sprintf('~/proj/att_awe/analysis/figures/fitCheck/%s',roi);
    print(fname,'-dpdf');
    if isfield(fit.fit,'scm'), fit.fit = rmfield(fit.fit,'scm'); end
    if isfield(fit.fit,'covar'), fit.fit = rmfield(fit.fit,'covar'); end

    fits{end+1} = fit;
end

%% Calculate Deconvolution by Condition

deconvo = {};
for ri = 1:length(rois)
    roi = rois{ri};       
    
    d = constructD(neural.tSeries.(name).(roi).tSeries,neural.SCM.(name).(roi).stimVol,0.5,20,neural.tSeries.(name).(roi).concatInfo,'none','deconv',0);
    decon = getr2timecourse(d.timecourse,d.nhdr,d.hdrlenTR,d.scm,d.framePeriod,d.verbose);
    decon = rmfield(decon,'scm');   
    decon = rmfield(decon,'covar');

    deconvo{end+1} = decon;
end

neural.(name).fits = fits;
neural.(name).deconvo = deconvo;















%%%%%%%%%%%%%%%%%%%%
%%   constructD   %%
%%%%%%%%%%%%%%%%%%%%
function d = constructD(timecourse,stimvol,framePeriod,hdrlen,concatInfo,option,fitType,verbose)

% make sure timecourse looks ok.
if (size(timecourse,2) == 1) && (size(timecourse,1) > 1)
  mrWarnDlg(sprintf('(fitTimecourse) Timecourse found to be %ix1 instead of 1x%i. Taking transpose',size(timecourse,1),size(timecourse,1)));
  timecourse = timecourse';
end

% first get stimulus convolution matrix for each condition separately
d.dim(4) = size(timecourse,2);
% note that this hdrlen is different from
% below, becuase it is in # of TR's not
% in seconds.
hdrlenTR = floor(hdrlen/framePeriod)+1;
d.hdrlen = hdrlenTR;
d.concatInfo = concatInfo;
d.nFrames = d.dim(4);

if verbose
  disppercent(-inf,'(fitTimecourse) Constructing stimulus convolution matrices for each condition');
end
for stimNum = 1:length(stimvol)
  d.stimvol{1} = stimvol{stimNum};
  d = makescm(d,[],0);
  eachSCM{stimNum} = d.scm;
  if verbose,disppercent(stimNum/length(stimvol));end
end
if verbose,disppercent(inf);end

% make an scm for the whole design
if verbose
  disppercent(-inf,'(fitTimecourse) Constructing stimulus convolution matrix for whole design');
end
d.stimvol = stimvol;
d = makescm(d,[],1);
allSCM = d.scm;
if verbose,disppercent(inf);end

% make an scm for each stimulus 
% now make a d structure to pass around
clear d;
d.dim = size(timecourse);
d.stimvol = stimvol;
d.framePeriod = framePeriod;
d.concatInfo = concatInfo;
d.eachSCM = eachSCM;
d.scm = allSCM;
d.hdrlen = hdrlen;
d.hdrlenTR = hdrlenTR;
d.nhdr = length(stimvol);
d.timecourse = timecourse;
d.applyFiltering = 1;
d.zeroMean = 0;
d.deconvModel = 1;
d.nFrames = max(d.dim);
d.verbose = verbose;