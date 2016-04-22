function [neural] = cc_performAnalysis(neural,name)

neural.fit = struct;

rois = neural.shortROIs;

if ~isdir(fullfile('~/proj/att_awe/analysis/figures/fitCheck/')), mkdir(fullfile('~/proj/att_awe/analysis/figures/fitCheck/')); end

%% Calculate Non-Lin Fits (for amplitude)

fitter = 'fitType=nonlin';
amper = 'amplitudeType=fit2';
fits = cell(1,length(rois));
for ri = 1:length(rois)
    roi = rois{ri};       
    
    if isfield(neural.tSeries.(name),roi) && isfield(neural.SCM.(name),roi)
        fit = fitTimecourse(neural.tSeries.(name).(roi).tSeries,neural.SCM.(name).(roi).stimVol,.5,'concatInfo',neural.tSeries.(name).(roi).concatInfo,fitter,amper);
        figure(1)
        fname = sprintf('~/proj/att_awe/analysis/figures/fitCheck/%s',roi);
        print(fname,'-dpdf');
        if isfield(fit.fit,'scm'), fit.fit = rmfield(fit.fit,'scm'); end
        if isfield(fit.fit,'covar'), fit.fit = rmfield(fit.fit,'covar'); end

        fits{ri} = fit;
    end
end

%% Calculate Deconvolution by Condition

deconvo = cell(1,length(rois));
for ri = 1:length(rois)
    roi = rois{ri};       
    
    if isfield(neural.tSeries.(name),roi) && isfield(neural.SCM.(name),roi)

        d = constructD(neural.tSeries.(name).(roi).tSeries,neural.SCM.(name).(roi).stimVol,0.5,20,neural.tSeries.(name).(roi).concatInfo,'none','deconv',0);
        decon = getr2timecourse(d.timecourse,d.nhdr,d.hdrlenTR,d.scm,d.framePeriod,d.verbose);
        decon = rmfield(decon,'scm');   
        decon = rmfield(decon,'covar');

        deconvo{ri} = decon;
    end
end

neural.(name).fits = fits;
neural.(name).deconvo = deconvo;














