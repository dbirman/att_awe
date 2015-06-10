function allData = ER_analysis(redo_all, redo_ts, redo_sv)

if redo_all
    redo_ts = 1;
    redo_sv = 1;
end

pre = '~/data/cohcon/s0300_pilot/s0300';
% folders = {'20150511'};
folders = {'20150509','20150513'};

thresh = .06;

curAnalysis = {'both_ER'};
ROIs = {{'l_v1','l_hmt','r_v1','r_hmt'}};
shortROIs = {{'v1','hmt'}};

total_analysis = sprintf('%s_%s.mat',strcat(folders{1}),strcat(curAnalysis{1}));

tempFolder = '~/data/temp';

tempSave = fullfile(tempFolder,total_analysis);

if isfile(tempSave)
    load(tempSave)
else
    data = struct;
    data.ROIs = ROIs;
    data.rROIs = shortROIs;
    data.analyses = curAnalysis;
    data.folders = folders;
    data.pre = pre;
end

if redo_ts || redo_sv
    for fi = 1:length(folders)
        folder = fullfile(sprintf('%s%s',pre,folders{fi}));
        fname = strcat('f',folders{fi});
        for ai = 1:length(curAnalysis)
            analysis = curAnalysis{ai};
            rois = ROIs{ai};

            mrQuit;

            if redo_all
                data.(fname) = struct;
                data.(fname).(analysis) = struct;
            end

            [tSeries,stimVol,concatInfo,stimNames] = loadERAnalysis(folder,rois,analysis, redo_ts, redo_sv, thresh);

            if redo_ts
                data.(fname).(analysis).tSeries = tSeries;
            end
            if redo_sv
                data.(fname).(analysis).stimVol = stimVol;
                data.(fname).(analysis).stimNames = stimNames;
            end
            data.(fname).(analysis).concatInfo = concatInfo;
        end
    end

    % At this point each cell is for each folder, and each sub-cell in tSeries
    % is for the ROIs. So we can run the analysis by concatenating the
    % appropriate runs, and then pulling out of the de-convolution the relevant
    % information.


    save(tempSave,'data');
end

data.rROIs = shortROIs;
allData = performAnalysis(data);

function allData = performAnalysis(data)

% Start by concatenating, within each ROI, within each analysis

% Note that the left/right hemisphere stimvols actually correspond to
% opposite mappings (i.e. lCoh or rCoh) of stimuli. But we can combine them
% by ignoring the prefix later on, so it's safe to concatenate everything
% here as long as it's coming form the same ROI.

for ai = 1:length(data.analyses)
    analysis = data.analyses{ai};
    rois = data.ROIs{ai};
    
    for ri = 1:length(rois)
        croi = rois{ri};
        % get the current ROI that we are in
        [~, rNum] = parseROI(croi);
        % get the ROI that this maps onto (we just ignore hemisphere)
        rroi = data.rROIs{ai}{rNum};
        
        if isfield(data,'concat') && isfield(data.concat,analysis) && isfield(data.concat.(analysis),rroi)
            % If the concat field already exists, we are going to add to it
            roiTSeries = data.concat.(analysis).(rroi).tSeries;
            roiStimVols = data.concat.(analysis).(rroi).stimvol;
            roiConcat = data.concat.(analysis).(rroi).concatInfo;
        else
            % If not, we're just going to start from scratch
            roiTSeries = [];
            roiStimVols = [];
            roiConcat = {};
        end
        
        % Now we do the actual concatentations
        for fi = 1:length(data.folders)
            folder = data.folders{fi};
            cdata = data.(strcat('f',folder)).(analysis);
            
            if isempty(roiTSeries)
                % If this is the first timeseries we just save the data
                roiTSeries = cdata.tSeries{ri};
                roiStimVols = cdata.stimVol{ri};
                roiConcat = cdata.concatInfo;
                % note stimnames is the same across all folders for an ROI
                % so this is okay. BUT: this only applies once prefixes
                % have been ignored obviously...
                data.concat.(analysis).(rroi).stimNames = cdata.stimNames{ri};
            else
                % Otherwise we concatenate to the previous data, we also
                % concatenate the stimVols even though the prefixes may be
                % different.
                [roiTSeries, roiStimVols, roiConcat] = concatRuns({roiTSeries, cdata.tSeries{ri}},{roiStimVols, cdata.stimVol{ri}},{roiConcat, cdata.concatInfo});
            end
        end
        % Save the data for the next run
        data.concat.(analysis).(rroi).tSeries = roiTSeries;
        data.concat.(analysis).(rroi).stimvol = roiStimVols;
        data.concat.(analysis).(rroi).concatInfo = roiConcat;
    end
end

allData = {};


cueds = {'uncued','cued'};
values = {'coherence','contrast'};

for ai = 1:length(data.analyses)
    analysis = data.analyses{ai};
    rois = data.rROIs{ai};
    for ri = 1:length(rois)
        roi = rois{ri};       
        [vals, cuedTask, valueType] = parseNames(data.concat.(analysis).(roi).stimNames);        
        fit = fitTimecourse(data.concat.(analysis).(roi).tSeries,data.concat.(analysis).(roi).stimvol,.5,'concatInfo',data.concat.(analysis).(roi).concatInfo,'fitType=deconv','amplitudeType=fit1');

%         fit = fitTimecourse(data.concat.(analysis).(roi).tSeries,data.concat.(analysis).(roi).stimvol,.5,'concatInfo',data.concat.(analysis).(roi).concatInfo,'fitType=glm','amplitudeType=area');
        
        
        for z = 1:length(fit.amplitude)
            if cuedTask(z) <= 2
                % ignore the 'catch' conditions for now
                amp = fit.amplitude(z); % amplitude
                ase = fit.amplitudeSTE(z);
                stim = valueType(z);
                val = vals(z);
                main = cuedTask(z);
                % since we know there are no catch trials, if cued==stim,
                % then it was attended
                cued = 1 + (stim == main);
                N = length(data.concat.(analysis).(roi).stimvol{z});
                
                if ~isfield(allData, roi), allData.(roi) = struct; end       
                if ~isfield(allData.(roi), values{stim}), allData.(roi).(values{stim}) = struct; end
                if ~isfield(allData.(roi).(values{stim}), cueds{cued})
                    allData.(roi).(values{stim}).(cueds{cued}) = struct;
                    allData.(roi).(values{stim}).(cueds{cued}).i = [];
                    allData.(roi).(values{stim}).(cueds{cued}).a = [];
                    allData.(roi).(values{stim}).(cueds{cued}).ase = [];
                    allData.(roi).(values{stim}).(cueds{cued}).N = [];
                    if isfield(fit,'canonicalResponse')
                        allData.(roi).(values{stim}).(cueds{cued}).canon = fit.canonicalResponse;
                    end
                    allData.(roi).(values{stim}).(cueds{cued}).times = 0.25:.5:20.25;
                end

                
                % intensity
                allData.(roi).(values{stim}).(cueds{cued}).i(end+1) = val;
                % amplitude
                allData.(roi).(values{stim}).(cueds{cued}).a(end+1) = amp;
                % error
                allData.(roi).(values{stim}).(cueds{cued}).ase(end+1) = ase;
                % N
                allData.(roi).(values{stim}).(cueds{cued}).N(end+1) = N;
            end
        end
    end
end

function [tSeries, stimVol, concatInfo,stimNames] = loadERAnalysis(folder,allROI,curA, r_ts, r_sv, thresh)
%%
cdir = pwd;
cd(folder);


%%
view = newView();
view = viewSet(view,'curGroup','Concatenation');
scans = viewGet(view,'nScans');
view = viewSet(view,'curScan',scans);
if r_ts
    view = loadAnalysis(view,sprintf('erAnal/%s',curA));
    analysis = viewGet(view,'analysis');
    d = analysis.d{1};
    d.scanNum = 1;
    d.groupNum = view.curGroup;
    d = loadroi(d,allROI);
end
concatInfo = viewGet(view,'concatInfo');

% generate stimVol cell
% format: {'rCon=XX,task=1', 'rCon=XX,task=2', etc...};'
stimVol = {}; stimNames = {};
sideVol{1} = [];sideVol{2} = []; sideNames{1} = [];sideNames{2} = [];
for ri = 1:length(allROI)
    roi = allROI{ri};
    if strfind(roi,'l_')
        prefix = 'r';
        side = 1;
    else
        prefix = 'l';
        side = 2;
    end

    if r_sv
        if isempty(sideVol{side})
            allStims = {};
            for t = 1:4
                for con = [.2 .4 .6 .8]
                    allStims{end+1} = {sprintf('%sCon=[%0.3f]',prefix,con), sprintf('nTask=[%i]',t)};
                end
                for coh = [0 .1 .25 .7]
                    allStims{end+1} = {sprintf('%sCoh=[%0.3f]',prefix,coh), sprintf('nTask=[%i]',t)};
                end
            end

            [stimVol{end+1}, stimNames{end+1}, ~] = getStimvol(view,allStims);
            sideVol{side} = stimVol{end};
            sideNames{side} = stimNames{end};
        else
            stimVol{end+1} = sideVol{side};
            stimNames{end+1} = stimNames{side};
        end
    else
        stimVol{end+1} = {};
        stimNames{end+1} = {};
    end
end


tSeries = {};
if r_ts
    scanDims = viewGet(view,'scanDims');
    r2 = viewGet(view,'overlayData',scans);
    for ri = 1:length(allROI)
        r = d.roi{ri};
        r.linearScanCoords = sub2ind(scanDims,r.scanCoords(1,:),r.scanCoords(2,:),r.scanCoords(3,:));

        r.r2 = r2(r.linearScanCoords);
        tSeries{end+1} = mean(r.tSeries(r.r2>thresh,:));
    end
end

clear view
clear analysis
clear d
clear r
clear r2

%%
cd(cdir);

function [side, roi] = parseROI(roiname)

sname = roiname(1:strfind(roiname,'_')-1);
rname = roiname(strfind(roiname,'_')+1:end);

if strcmp(sname,'l')
    side = 1;
else
    side = 2;
end

if strcmp(rname,'v1')
    roi = 1;
else
    roi = 2;
end


function [val, cuedTask, valueType] = parseNames(stimNames)
%% Note that parse names ignores the prefix entirely, so it's hemisphere independent
val = []; cuedTask = []; valueType = [];

vstr = 'Co';

for i = 1:length(stimNames)
    name = stimNames{i};
    % t can be 1, 2, (coherence, contrast) or 3, 4 (catch coherence, catch
    % contrast)
    t = str2num(name(strfind(name,'ask=')+4:end-1));
    % v is the value of coh/con that was shown on the screen on this trial
    v = str2num(name(strfind(name,vstr)+4:strfind(name,' &')));
    val(end+1) = v;
    cuedTask(end+1) = t;
    if strfind(name,'Con')
        valueType(end+1) = 2;
    else
        valueType(end+1) = 1;
    end
end