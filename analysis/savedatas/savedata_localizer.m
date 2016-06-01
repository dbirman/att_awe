function savedata_localizer(cfolder)

% Generic function to save data from localizer runs

% (1) Extract the left and right flat map voxel timecourses from the
% concatenation timeseries
% (2) Generate the design by loading the contrast, coherence, and timing
% responses, using the basecon/basecoh of the run (default timing=5)
% This is the ROI_data structure and gets saved in ~/data/cohcon_localizer

%% Move to Data WD
mrQuit
cd(fullfile('~/data/cohcon_localizer/',cfolder));

%% Testing
[session, groups, ~] = loadSession(pwd);
%% Setup a view + Load Concatenation
view = newView();
view = viewSet(view,'curGroup','Concatenation');
view = viewSet(view,'curScan',1); % make sure scan # is correct

%% Get the sort order from the _all_ analysis
view = loadAnalysis(view,sprintf('erAnal/%s','all')); % check analysis name!
r2 = viewGet(view,'overlayData',view.curScan);

%% Pull Stimvols
[convol, conconds, ~] = getStimvol(view,'contrast','taskNum=1','phaseNum=2');
[cohvol, cohconds, ~] = getStimvol(view,'coherence','taskNum=1','phaseNum=2');
[timvol, timconds, ~] = getStimvol(view,'timing','taskNum=1','phaseNum=2');
[conval,cohval,timval] = parseIndivNames(conconds,cohconds,timconds);
timval = timval*2; % use the .5s versions
timdef = 5; % default of 2.5s stim
basecon = 0.25;
basecoh = 0;

%% Convert to Long Form
condata = volval2long(convol,conval);
cohdata = volval2long(cohvol,cohval);
timdata = volval2long(timvol,timval);

%% Load occipital ROIs
rois = loadROITSeries(view,{'V1'},view.curScan,view.curGroup,'keepNAN=true');
rois = getSortIndex(view,rois,r2);

%% Generate design
% stimvol - basecon - newcon - basecoh - newcoh - timing
design = zeros(10000,6);
count = 1;

for vol = 1:size(condata,1)
    % take this contrast trial and find its coh/tim
    csv = condata(vol,1);
    ccon = condata(vol,2);
    cohidx = find(cohdata(:,1)==csv);
    if ~isempty(cohidx)
        ccoh = cohdata(cohidx,2);
    else
        ccoh = basecoh;
    end
    timidx = find(timdata(:,1)==csv);
    if ~isempty(timidx)
        ctim = timdata(timidx,2);
    else
        ctim = timdef;
    end
    design(count,:) = [csv basecon ccon basecoh ccoh ctim];
    count = count + 1;
end
for vol = 1:size(cohdata,1)
    csv = cohdata(vol,1);
    if ~any(design(:,1)==csv)
        ccoh = cohdata(vol,2);
        % we didn't put this trial in, check for timing
        timidx = find(timdata(:,1)==csv);
        if ~isempty(timidx)
            ctim = timdata(timidx,2);
        else
            ctim = timdef;
        end
        design(count,:) = [csv basecon basecon basecoh ccoh ctim];
        count = count+1;
    end
end
for vol = 1:size(timdata,1)
    csv = timdata(vol,1);
    if ~any(design(:,1)==csv)
        warning('failure');
    end
end
design = design(1:count-1,:);

%% Save data
savefolder = '~/data/cohcon_localizer';
fname = sprintf('%s_data.mat',cfolder);
fname = fullfile(savefolder,fname);

data.rois = rois;
data.design = design;
data.basecon = basecon;
data.basecoh = basecoh;
data.concatInfo = viewGet(view,'concatInfo');
data.TR = 0.5;

save(fname,'data');

clear data
