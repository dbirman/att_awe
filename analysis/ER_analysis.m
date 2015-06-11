function allData = ER_analysis(redo_all, redo_ts, redo_sv)

if redo_all
    redo_ts = 1;
    redo_sv = 1;
end

pre = '~/data/cohcon/s0300_pilot/s0300';
% folders = {'20150509'};
folders = {'20150509','20150511','20150513'};

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
    f = figure;
    for ri = 1:length(rois)
        roi = rois{ri};       
        [conVals, cohVals, cuedTask] = parseNames(data.concat.(analysis).(roi).stimNames); 
        fitter = 'fitType=glm';
        amper = 'amplitudeType=area';
        fit = fitTimecourse(data.concat.(analysis).(roi).tSeries,data.concat.(analysis).(roi).stimvol,.5,'concatInfo',data.concat.(analysis).(roi).concatInfo,fitter,amper);

        figure(f)
        N = [];
        
        %% Do matrix calculation
        for si = 1:length(data.concat.(analysis).(roi).stimvol)
            N(end+1) = length(data.concat.(analysis).(roi).stimvol{si});
        end
%         fit = fitTimecourse(data.concat.(analysis).(roi).tSeries,data.concat.(analysis).(roi).stimvol,.5,'concatInfo',data.concat.(analysis).(roi).concatInfo,'fitType=glm','amplitudeType=area');
        conVals = conVals(cuedTask<=2);
        cohVals = cohVals(cuedTask<=2);
        amps = fit.amplitude(cuedTask<=2);
        ampse = fit.amplitudeSTE(cuedTask<=2);
        N = N(cuedTask<=2);
        cuedTask = cuedTask(cuedTask<=2);
        
        ucon = unique(conVals);
        ucoh = unique(cohVals);
        
        mat{1} = zeros(length(ucon),length(ucoh));
        mat{2} = zeros(length(ucon),length(ucoh));
        matse{1} = zeros(length(ucon),length(ucoh));
        matse{2} = zeros(length(ucon),length(ucoh));
        
        for ai = 1:length(amps)
            mat{cuedTask(ai)}(find(ucon==conVals(ai)),find(ucoh==cohVals(ai))) = amps(ai);
            matse{cuedTask(ai)}(find(ucon==conVals(ai)),find(ucoh==cohVals(ai))) = ampse(ai);
        end
        
        
        con = [.2 .4 .6 .8];
        coh = [0 .1 .25 .7];
        %% X axis = contrast
        subplot(length(rois),2,(ri-1)*2+1); hold on
        clist1 = brewermap(4,'Oranges');
        clist2 = brewermap(4,'Greens');
        for i = 1:4
            plot(con,mat{1}(i,:),'-','Color',clist1(i,:));
            errorbar(con,mat{1}(i,:),matse{1}(i,:),'*','Color',clist1(i,:));
            plot(con,mat{2}(i,:),'-','Color',clist2(i,:));
            errorbar(con,mat{2}(i,:),matse{2}(i,:),'*','Color',clist2(i,:));
        end
        title(sprintf('%s: Orange = Cued Coh, Green = Cued Con.',roi));
        xlabel('Contrast');
        ylabel('Response Amplitude');
        axis([0 1 0 450])
        
        %% X axis = coherence
        subplot(length(rois),2,(ri-1)*2+2); hold on
        clist1 = brewermap(4,'Oranges');
        clist2 = brewermap(4,'Greens');
        for i = 1:4
            m1 = mat{1}';
            m2 = mat{2}';
            plot(coh,m1(i,:),'-','Color',clist1(i,:));
            errorbar(coh,m1(i,:),matse{1}(i,:),'*','Color',clist1(i,:));
            plot(coh,m2(i,:),'-','Color',clist2(i,:));
            errorbar(coh,m2(i,:),matse{2}(i,:),'*','Color',clist2(i,:));
        end
        title(sprintf('%s: Orange = Cued Coh, Green = Cued Con.',roi));
        xlabel('Coherence');
        ylabel('Response Amplitude');
        axis([0 1 0 450])
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
view = viewSet(view,'curScan',1);
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
            for t = 1:2
                for con = [.2 .4 .6 .8]
                    for coh = [0 .1 .25 .7]
                        allStims{end+1} = {sprintf('%sCon=[%0.3f]',prefix,con),sprintf('%sCoh=[%0.3f]',prefix,coh), sprintf('nTask=[%i]',t)};
                    end
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


function [conVal, cohVal, cuedTask] = parseNames(stimNames)
%% Note that parse names ignores the prefix entirely, so it's hemisphere independent
conVal = []; cohVal = []; cuedTask = [];

for i = 1:length(stimNames)
    name = stimNames{i};
    ands = strfind(name,' &');
    
    if strfind(name,'Con=')
        conVal(end+1) = str2num(name(strfind(name,'Con=')+4:ands(1))); % note contrast first
    end
    if strfind(name,'Coh=')
        cohVal(end+1) = str2num(name(strfind(name,'Coh=')+4:ands(2)-2));
    end
    % t can be 1, 2, (coherence, contrast) or 3, 4 (catch coherence, catch
    % contrast)
    cuedTask(end+1) = str2num(name(strfind(name,'ask=')+4:end-1));
end