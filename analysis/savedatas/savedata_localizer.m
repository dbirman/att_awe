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

%% Setup a view + Load Concatenation
view = newView();
view = viewSet(view,'curGroup','Concatenation');
view = viewSet(view,'curScan',1); % make sure scan # is correct

%% Get the mean timeseries using the reversed pRF
view = loadAnalysis(view,sprintf('erAnal/%s','all')); % check analysis name!
ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','LO1','LO2','MT'};
% ROIs = {'V1'};
pfxs = {'l','r'};
allROIs = {};
for ri = 1:length(ROIs)
    for pi = 1:length(pfxs)
        allROIs{end+1} = strcat(pfxs{pi},ROIs{ri});
    end
end
analysis = view.analyses{1};
rois = loadROITSeries(view,allROIs,view.curScan,view.curGroup,'keepNAN=true');

tSeries = cell(1,length(rois)); % meaned tSeries with prf
rtSeries = cell(1,length(rois));
r2s = cell(1,length(rois));
rights = cell(1,length(rois));
lefts = cell(1,length(rois));
corrs = cell(1,length(rois));

scanDims = viewGet(view,'scanDims');

left = []; right = []; r2 = [];
% check which overlay is which to make sure they get sorted properly
for i = 1:length(analysis.overlays)
    cOverlay = analysis.overlays(i);
    if strfind(cOverlay.name,'r2')
        r2 = cOverlay.data{1};
        disp(sprintf('(sd_loc) Setting overlay %i to r2',i));
    elseif strfind(cOverlay.name,'left')
        left = cOverlay.data{1};
        disp(sprintf('(sd_loc) Setting overlay %i to left',i));
    elseif strfind(cOverlay.name,'right')
        right = cOverlay.data{1};
        disp(sprintf('(sd_loc) Setting overlay %i to right',i));
    else
        warning('failure');
        keyboard
    end
end

if isempty(left) || isempty(right)
    disp('(sd_loc) Failed to find the overlays, averaging across entire ROI');
end

% just incase we want this at some point?
rois = getSortIndex(view,rois,r2);

for ri = 1:length(rois)
    r = rois{ri};
    r.linearScanCoords = sub2ind(scanDims,r.scanCoords(1,:),r.scanCoords(2,:),r.scanCoords(3,:));
    
    cr2 = r2(r.linearScanCoords);   
    if ~isempty(left) && ~isempty(right) 
        rightO = right(r.linearScanCoords);
        leftO = left(r.linearScanCoords);
        rightO(isnan(rightO))=0;
        leftO(isnan(leftO))=0;

        idxs = ~any(isnan(r.tSeries),2);
        if any(~idxs)
            warning('Failure');
            keyboard
        end
    %     tSeriesnoNaN = r.tSeries;
    %     tSeriesnoNaN(isnan(tSeriesnoNaN)) = 0;

        r2s{ri} = cr2;
        rights{ri} = rightO;
        lefts{ri} = leftO;

        if strcmp(r.name(1),'l')
            tSeries{ri} = (rightO*r.tSeries)/sum(rightO);
            corrs{ri} = corr(cr2',rightO');
        else
            tSeries{ri} = (leftO*r.tSeries)/sum(leftO);
            corrs{ri} = corr(cr2',leftO');
        end
    else
        tSeries{ri} = mean(r.tSeries);
        corrs{ri} = 0;
        rights{ri} = [];
        lefts{ri} = [];
        r2s{ri} = cr2;
    end
    r2cutoff = 0.5;
    while sum(cr2>r2cutoff)<25
        r2cutoff = r2cutoff-0.001;
    end
    rtSeries{ri} = mean(r.tSeries(cr2>r2cutoff,:));
end

%% testing
% r = rois{19};    r.linearScanCoords = sub2ind(scanDims,r.scanCoords(1,:),r.scanCoords(2,:),r.scanCoords(3,:));
% cr2 = r2(r.linearScanCoords);
% clf, hold on
% r2cutoff = 0.5;
% while sum(cr2>r2cutoff)<25
%     r2cutoff = r2cutoff-0.001;
% end
% point5 = mean(r.tSeries(cr2>r2cutoff,:));
% rightO = left(r.linearScanCoords);
% rightO = rightO;
% rightO(isnan(rightO))=0;
% rightf = (rightO*r.tSeries)/sum(rightO);
% plot(point5(1000:2200));
% plot(rightf(1000:2200),'r');
% a = axis;
% axis([a(1) a(2) .95 1.05]);

%% Pull Stimvols
[lConvol, lConconds, ~] = getStimvol(view,'lCon','taskNum=1','phaseNum=2');
[rConvol, rConconds, ~] = getStimvol(view,'rCon','taskNum=1','phaseNum=2');
[lCohvol, lCohconds, ~] = getStimvol(view,'lCoh','taskNum=1','phaseNum=2');
[rCohvol, rCohconds, ~] = getStimvol(view,'rCoh','taskNum=1','phaseNum=2');
[convol, conconds, ~] = getStimvol(view,'contrast','taskNum=1','phaseNum=2');
[cohvol, cohconds, ~] = getStimvol(view,'coherence','taskNum=1','phaseNum=2');
[timvol, timconds, ~] = getStimvol(view,'timing','taskNum=1','phaseNum=2');
[taskvol, taskconds, ~] = getStimvol(view,'task','taskNum=1','phaseNum=2');
[correctvol, correctconds, ~] = getStimvol(view,'correct','taskNum=1','phaseNum=2');
if ~isempty(correctvol)
    warning('ADD CORRECT CHECKS!');
    stop = 1;
end
[conval,cohval,timval,lconval,rconval,lcohval,rcohval,taskval] = parseIndivNames(conconds,cohconds,timconds,lConconds,rConconds,lCohconds,rCohconds,taskconds);
timval = timval*2; % use the .5s versions
basecon = 0.25;
basecoh = 0;

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
end

% we won't use condata/cohdata again, we will only explicitly model
% lCon/rCon and lCoh/rCoh since the ROIs are by hemi anyways

%% Grab all the unique stimvols so we know how many trials we should get
allsv = unique([lcondata(:,1);rcondata(:,1);rcohdata(:,1);lcohdata(:,1);timdata(:,1);taskdata(:,1)]);

%% Generate design
design = zeros(10000,9);
count = 1;

for ci = 1:length(allsv)
    csv = allsv(ci);
    % this lCon event, find its rCon, lCoh, rCoh, tim, and task
    conidx = condata(:,1)==csv;
    cohidx = cohdata(:,1)==csv;
    lConidx = find(lcondata(:,1)==csv);
    if isempty(lConidx)
        lCon = condata(conidx,2);
    else
        lCon = lcondata(lConidx,2);
    end
    % get rCon
    rConidx = find(rcondata(:,1)==csv);
    if isempty(rConidx)
        rCon = condata(conidx,2);
    else
        rCon = rcondata(rConidx,2);
    end
    % get lCoh
    lCohidx = find(lcohdata(:,1)==csv);
    if isempty(lCohidx)
        lCoh = cohdata(cohidx,2);
    else
        lCoh = lcohdata(lCohidx,2);
    end
    % get rCoh
    rCohidx = find(rcohdata(:,1)==csv);
    if isempty(rCohidx)
        rCoh = cohdata(cohidx,2);
    else
        rCoh = rcohdata(rCohidx,2);
    end
    % get tim
    timidx = find(timdata(:,1)==csv);
    if ~isempty(timidx)
        tim = timdata(timidx,2);
    else
        % figure out what kind of trial this is (cohxcon or cohxcon+att)
        if lCon==rCon
            tim = 5;
        else
            tim = 1;
        end
    end
    % get task
    taskidx = find(taskdata(:,1)==csv);
    if ~isempty(taskidx)
        task = taskdata(taskidx,2);
    else
        % no task information, so this trial must be a fixation task
        if lCon==rCon
            task = 0;
        else
            warning('Shouldn''t be able to get here...');
            keyboard
        end
    end
    design(count,:) = [csv basecon lCon rCon basecoh lCoh rCoh tim task];
    % stimvol basecon lcon rcon basecoh lcoh rcoh timing task
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
data.rtSeries = rtSeries;
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
data.concatInfo = concatInfo;
data.TR = 0.5;

save(fname,'data');

clear data

disp(sprintf('Saving data for %s completed successfully',cfolder));