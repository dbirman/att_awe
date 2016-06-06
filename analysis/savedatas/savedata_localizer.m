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
    
    idxs = ~any(isnan(r.tSeries),2);
    if any(~idxs)
        keyboard
    end
%     tSeriesnoNaN = r.tSeries;
%     tSeriesnoNaN(isnan(tSeriesnoNaN)) = 0;
    
    rightMean = (rightO*r.tSeries)/sum(rightO);
    leftMean = (leftO*r.tSeries)/sum(leftO);
    
    r2s{ri} = cr2;
    rights{ri} = rightO;
    lefts{ri} = leftO;
    
    if strcmp(r.name(1),'l')
        tSeries{ri} = rightMean;
    else
        tSeries{ri} = leftMean;
    end
end

%% Pull Stimvols
[lConvol, lConconds, ~] = getStimvol(view,'lCon','taskNum=1','phaseNum=2');
[rConvol, rConconds, ~] = getStimvol(view,'rCon','taskNum=1','phaseNum=2');
[convol, conconds, ~] = getStimvol(view,'contrast','taskNum=1','phaseNum=2');
[cohvol, cohconds, ~] = getStimvol(view,'coherence','taskNum=1','phaseNum=2');
[timvol, timconds, ~] = getStimvol(view,'timing','taskNum=1','phaseNum=2');
[conval,cohval,timval] = parseIndivNames(conconds,cohconds,timconds);
timval = timval*2; % use the .5s versions
timdef = 5; % default of 2.5s stim
basecon = 0.25;
basecoh = 0;

if ~isempty(lConvol)
    warning('Write the code for lCoh/rCoh and Con');
    keyboard
end

%% Convert to Long Form
condata = volval2long(convol,conval);
cohdata = volval2long(cohvol,cohval);
timdata = volval2long(timvol,timval);

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

data.tSeries = tSeries;
data.ROIs = allROIs;
data.roi.r2 = r2s;
data.roi.right = rights;
data.roi.left = lefts;
data.design = design;
data.basecon = basecon;
data.basecoh = basecoh;
data.concatInfo = viewGet(view,'concatInfo');
data.TR = 0.5;

save(fname,'data');

clear data
