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

%% Get the mean timeseries using the reversed pRF
view = loadAnalysis(view,sprintf('erAnal/%s','all')); % check analysis name!
ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','LO1','LO2','MT'};
pfxs = {'l','r'};
allROIs = {};
for ri = 1:length(ROIs)
    for pi = 1:length(pfxs)
        allROIs{end+1} = strcat(pfxs{pi},ROIs{ri});
    end
end
analysis = view.analyses{1};
rois = loadROITSeries(view,allROIs,view.curScan,view.curGroup,'keepNAN=true');

tSeries = cell(1,length(rois)); % meaned tSeries
r2s = cell(1,length(rois));
rights = cell(1,length(rois));
lefts = cell(1,length(rois));
corrs = cell(1,length(rois));

scanDims = viewGet(view,'scanDims');

warning('Should check which overlay is which here, not just assume order');
r2 = analysis.overlays(1).data{1};
rois = getSortIndex(view,rois,r2);
right = analysis.overlays(2).data{1};
left = analysis.overlays(3).data{1};
for ri = 1:length(rois)
    r = rois{ri};
    r.linearScanCoords = sub2ind(scanDims,r.scanCoords(1,:),r.scanCoords(2,:),r.scanCoords(3,:));
    
    cr2 = r2(r.linearScanCoords);    
    rightO = right(r.linearScanCoords);
    leftO = left(r.linearScanCoords);
    rightO(isnan(rightO))=0;
    leftO(isnan(leftO))=0;
    rightO = rightO.^2;
    leftO = leftO.^2;
    
    idxs = ~any(isnan(r.tSeries),2);
    if any(~idxs)
        keyboard
    end
%     tSeriesnoNaN = r.tSeries;
%     tSeriesnoNaN(isnan(tSeriesnoNaN)) = 0;
    
    r2s{ri} = cr2;
    rights{ri} = rightO;
    lefts{ri} = leftO;
    
    if strcmp(r.name(1),'l')
        tSeries{ri} = (rightO*r.tSeries)/sum(rightO);
        corrs{ri} = corr(cr2,rightO);
    else
        tSeries{ri} = (leftO*r.tSeries)/sum(leftO);
        corrs{ri} = corr(cr2,leftO);
    end
end

%% Pull Stimvols
[lConvol, lConconds, ~] = getStimvol(view,'lCon','taskNum=1','phaseNum=2');
[rConvol, rConconds, ~] = getStimvol(view,'rCon','taskNum=1','phaseNum=2');
[lCohvol, lCohconds, ~] = getStimvol(view,'lCoh','taskNum=1','phaseNum=2');
[rCohvol, rCohconds, ~] = getStimvol(view,'rCoh','taskNum=1','phaseNum=2');
[convol, conconds, ~] = getStimvol(view,'contrast','taskNum=1','phaseNum=2');
[cohvol, cohconds, ~] = getStimvol(view,'coherence','taskNum=1','phaseNum=2');
[timvol, timconds, ~] = getStimvol(view,'timing','taskNum=1','phaseNum=2');
[taskvol, taskconds, ~] = getStimvol(view,'task','taskNum=1','phaseNum=2');
[conval,cohval,timval,lconval,rconval,lcohval,rcohval,taskval] = parseIndivNames(conconds,cohconds,timconds,lConconds,rConconds,lCohconds,rCohconds,taskconds);
timval = timval*2; % use the .5s versions
basecon = 0.25;
basecoh = 0;
taskdef = -1;

%% Convert to Long Form
condata = volval2long(convol,conval);
cohdata = volval2long(cohvol,cohval);
timdata = volval2long(timvol,timval);
lcondata = volval2long(lConvol,lconval);
rcondata = volval2long(rConvol,rconval);
lcohdata = volval2long(lCohvol,lcohval);
rcohdata = volval2long(rCohvol,rcohval);
taskdata = volval2long(taskvol,taskval);

if isempty(lcondata)
    % if lcon is empty all four are probably empty
    lcondata = condata;
    rcondata = condata;
    lcohdata = cohdata;
    rcohdata = cohdata;
    % default to the 2.5 s stimulus and fixation task
    timdef = 5;
    taskdef = 0;
else
    % if we have lCon/rCon default to the .5s stimulus
    timdef = 1;
end

% we won't use condata/cohdata again, we will only explicitly model
% lCon/rCon and lCoh/rCoh since the ROIs are by hemi anyways

%% Grab all the unique stimvols so we know how many trials we should get
allsv = unique([lcondata(:,1);rcondata(:,1);rcohdata(:,1);lcohdata(:,1);timdata(:,1);taskdata(:,1)]);

%% Generate design
% stimvol - basecon - newcon - basecoh - newcoh - timing - attention 0/1/2
%                                               0=fixation, 1=coh, 2=con
design = zeros(10000,9);
count = 1;

for vol = 1:size(lcondata,1)
    % this lCon event, find its rCon, lCoh, rCoh, tim, and task
    csv = lcondata(vol,1);
    lCon = lcondata(vol,2);
    % get rCon
    rConidx = find(rcondata(:,1)==csv);
    if isempty(rConidx), warning('Failure'); keyboard; end
    rCon = rcondata(rConidx,2);
    % get lCoh
    lCohidx = find(lcohdata(:,1)==csv);
    if isempty(lCohidx), warning('Failure'); keyboard; end
    lCoh = lcohdata(lCohidx,2);
    % get rCoh
    rCohidx = find(rcohdata(:,1)==csv);
    if isempty(rCohidx), warning('Failure'); keyboard; end
    rCoh = rcohdata(rCohidx,2);
    % get tim
    timidx = find(timdata(:,1)==csv);
    if ~isempty(timidx)
        tim = timdata(timidx,2);
    else
        tim = timdef;
    end
    % get task
    taskidx = find(taskdata(:,1)==csv);
    if ~isempty(taskidx)
        task = taskdata(taskidx,2);
    else
        if taskdef==-1, warning('Failure'); keyboard; end
        task = taskdef;
    end
    design(count,:) = [csv basecon lCon rCon basecoh lCoh rCoh tim task];
    count = count+1;
end
design = design(1:count-1,:);

if size(design,1)~=size(allsv,1)
    warning('Some trials are missing, or too many trials were found');
    keyboard
end

%% Save data
savefolder = '~/data/cohcon_localizer';
fname = sprintf('%s_data.mat',cfolder);
fname = fullfile(savefolder,fname);

data.tSeries = tSeries;
data.ROIs = allROIs;
data.roi.r2 = r2s;
data.roi.right = rights;
data.roi.left = lefts;
data.roi.corrs = corrs;
data.design = design;
data.basecon = basecon;
data.basecoh = basecoh;
concatInfo = viewGet(view,'concatInfo');
data.runtrans = concatInfo.runTransition;
data.TR = 0.5;

save(fname,'data');

clear data

disp(sprintf('Saving data for %s completed successfully',cfolder));