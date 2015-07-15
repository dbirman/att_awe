function allData = ER_analysis(redo_all, redo_ts, redo_sv)

if redo_all
    redo_ts = 1;
    redo_sv = 1;
end

pre = '~/data/cohcon/s0300_pilot/s0300';
% folders = {'20150509'};
folders = {'20150509','20150511','20150513','20150618'};

thresh = .07;

curAnalysis = {'both_ER'};
ROIs = {{'l_v1','l_v2v','l_v2d','l_v3v','l_v3d','l_hmt','r_v1','r_v2v','r_v2d','r_v3v','r_v3d','r_hmt'}};
shortROIs = {{'v1','v2','v3','hmt'}};

total_analysis = sprintf('%s_%s.mat',strcat(folders{1}),strcat(curAnalysis{1}));

tempFolder = '~/data/temp';

tempSave = fullfile(tempFolder,total_analysis);

if isfile(tempSave) && ~redo_all
    load(tempSave)
else
    data = struct;
end

data.ROIs = ROIs;
data.rROIs = shortROIs;
data.analyses = curAnalysis;
data.folders = folders;
data.pre = pre;

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
%%
for ai = 1:length(data.analyses)
    analysis = data.analyses{ai};
    rois = data.ROIs{ai};
    
    for ri = 1:length(rois)
        croi = rois{ri};
        % get the current ROI that we are in
        [side, rNum] = parseROI(croi);
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
                if side==1
                    roiStimVols = cdata.stimVol{1};
                    data.concat.(analysis).(rroi).stimNames = cdata.stimNames{1};
                else
                    roiStimVols = cdata.stimVol{end};
                    data.concat.(analysis).(rroi).stimNames = cdata.stimNames{end};
                end
                roiConcat = cdata.concatInfo;
                % note stimnames is the same across all folders for an ROI
                % so this is okay. BUT: this only applies once prefixes
                % have been ignored obviously...
            else
                % Otherwise we concatenate to the previous data, we also
                % concatenate the stimVols even though the prefixes may be
                % different.
                if side==1
                    [roiTSeries, roiStimVols, roiConcat] = concatRuns({roiTSeries, cdata.tSeries{ri}},{roiStimVols, cdata.stimVol{1}},{roiConcat, cdata.concatInfo});
                else
                    [roiTSeries, roiStimVols, roiConcat] = concatRuns({roiTSeries, cdata.tSeries{ri}},{roiStimVols, cdata.stimVol{end}},{roiConcat, cdata.concatInfo});
                end
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

    %% Calculate Fits
for ai = 1:length(data.analyses)
    
    analysis = data.analyses{ai};
    rois = data.rROIs{ai};
    fits = {};
    for ri = 1:length(rois)
        roi = rois{ri};       
        [conVals, cohVals, cuedTask] = parseNames(data.concat.(analysis).(roi).stimNames); 
        fitter = 'fitType=nonlin';
        amper = 'amplitudeType=fit2';
        fits{end+1} = fitTimecourse(data.concat.(analysis).(roi).tSeries,data.concat.(analysis).(roi).stimvol,.5,'concatInfo',data.concat.(analysis).(roi).concatInfo,fitter,amper);
    end
    
        
    %% Figures
    f1 = figure;
    f2 = figure;
    for ri = 1:length(rois)
        roi = rois{ri};     
        fit = fits{ri};
        N = [];
        for si = 1:length(data.concat.(analysis).(roi).stimvol)
            N(end+1) = length(data.concat.(analysis).(roi).stimvol{si});
        end
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
        Nm{1} = zeros(length(ucon),length(ucoh));
        Nm{2} = zeros(length(ucon),length(ucoh));
        
        for ampi = 1:length(amps)
            mat{cuedTask(ampi)}(find(ucon==conVals(ampi)),find(ucoh==cohVals(ampi))) = amps(ampi);
            matse{cuedTask(ampi)}(find(ucon==conVals(ampi)),find(ucoh==cohVals(ampi))) = ampse(ampi);
            Nm{cuedTask(ampi)}(find(ucon==conVals(ampi)),find(ucoh==cohVals(ampi))) = N(ampi);
        end
        
        figure(f1)
        % X axis = contrast
        subplot(length(rois),2,(ri-1)*2+1); hold on
        clist1 = brewermap(4,'Oranges');
        clist2 = brewermap(4,'Purples');
        legVals = {};
        for i = 1:4
            plot(ucon,mat{1}(:,i),'-','Color',clist1(i,:));
            legVals{end+1} = sprintf('Cued Coh, Coh = %0.2f',ucoh(i));
        end
        for i = 1:4
            plot(ucon,mat{2}(:,i),'-','Color',clist2(i,:));
            legVals{end+1} = sprintf('Cued Con, Coh = %0.2f',ucoh(i));
        end
        legend(legVals);
        for i = 1:4
            errorbar(ucon,mat{1}(:,i),matse{1}(:,i),'*','Color',clist1(i,:));
            errorbar(ucon,mat{2}(:,i),matse{2}(:,i),'*','Color',clist2(i,:));
        end
        title(sprintf('%s: Orange = Cued Coh, Green = Cued Con.',roi));
        xlabel('Contrast');
        ylabel('Response Amplitude');
        axis([0 1 1 3])
        
        % X axis = coherence
        subplot(length(rois),2,(ri-1)*2+2); hold on
        for i = 1:4
            plot(ucoh,mat{1}(i,:),'-','Color',clist1(i,:));
            errorbar(ucoh,mat{1}(i,:),matse{1}(i,:),'*','Color',clist1(i,:));
            plot(ucoh,mat{2}(i,:),'-','Color',clist2(i,:));
            errorbar(ucoh,mat{2}(i,:),matse{2}(i,:),'*','Color',clist2(i,:));
        end
        title(sprintf('%s: Orange = Cued Coh, Green = Cued Con.',roi));
        xlabel('Coherence');
        ylabel('Response Amplitude');
        axis([0 1 1 3])
        
        
        figure(f2)
        subplot(length(rois),2,(ri-1)*2+1); hold on
        plot(ucon,mean(mat{1},2),'-','Color',clist1(end,:));
        plot(ucon,mean(mat{2},2),'-','Color',clist2(end,:));
        legend({'Cued Coh','Cued Con'});
        % get the pooled se...
        pooledse{1} = sum(matse{1}.*(Nm{1}-1),2) ./ sum(Nm{1}-1,2);
        pooledse{2} = sum(matse{2}.*(Nm{2}-1),2) ./ sum(Nm{1}-1,2);
        % pooled
        errorbar(ucon,mean(mat{1},2),1.96*pooledse{1},'*','Color',clist1(end,:));
        errorbar(ucon,mean(mat{2},2),1.96*pooledse{2},'*','Color',clist2(end,:));
        
        title(sprintf('%s: Orange = Cued Coh, Green = Cued Con.',roi));
        xlabel('Contrast');
        ylabel('Response Amplitude');
        axis([0 1 1 3])
        subplot(length(rois),2,(ri-1)*2+2); hold on
        plot(ucoh,mean(mat{1},1),'-','Color',clist1(end,:));
        plot(ucoh,mean(mat{2},1),'-','Color',clist2(end,:));
        legend({'Cued Coh','Cued Con'});
        % get the pooled se...
        pooledse{1} = sum(matse{1}.*(Nm{1}-1),1) ./ sum(Nm{1}-1,1);
        pooledse{2} = sum(matse{2}.*(Nm{2}-1),1) ./ sum(Nm{1}-1,1);
        % pooled
        errorbar(ucoh,mean(mat{1},1),1.96*pooledse{1},'*','Color',clist1(end,:));
        errorbar(ucoh,mean(mat{2},1),1.96*pooledse{2},'*','Color',clist2(end,:));
        title(sprintf('%s: Orange = Cued Coh, Green = Cued Con.',roi));
        xlabel('Coherence');
        ylabel('Response Amplitude');
        axis([0 1 1 3])
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
view = viewSet(view,'curScan',2);
if r_ts
    view = loadAnalysis(view,sprintf('erAnal/%s',curA));
    analysis = viewGet(view,'analysis');
    d = analysis.d{2};
    d.scanNum = 2;
    d.groupNum = view.curGroup;
    d = loadroi(d,allROI);
end
concatInfo = viewGet(view,'concatInfo');

% generate stimVol cell
% format: {'rCon=XX,task=1', 'rCon=XX,task=2', etc...};'
stimVol = {{},{}}; stimNames = {{},{}};
for ri = 1:length(allROI)
    roi = allROI{ri};
    if strfind(roi,'l_')
        prefix = 'r';
    else
        prefix = 'l';
    end

    if r_sv
        if strfind(roi,'l_')
            side = 1;
        else
            side = 2;
        end
        if ~isempty(allStims{side})
            allStims = {};
            for t = 1:2
                for con = [.2 .4 .6 .8]
                    for coh = [0 .1 .25 .7]
                        allStims{end+1} = {sprintf('%sCon=[%0.3f]',prefix,con),sprintf('%sCoh=[%0.3f]',prefix,coh), sprintf('nTask=[%i]',t)};
                    end
                end
            end

            [stimVol{side}, stimNames{side}, ~] = getStimvol(view,allStims);
        end
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

%ROIs = {{'l_v1','l_v2v','l_v2d','l_v3v','l_v3d','l_hmt','r_v1','r_v2v','r_v2d','r_v3v','r_v3d','r_hmt'}};


switch rname
    case 'v1'
        roi = 1;
    case 'v2v'
        roi = 2;
    case 'v2d'
        roi = 2;
    case 'v3v'
        roi = 3;
    case 'v3d'
        roi = 3;
    case 'hmt'
        roi = 4;
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