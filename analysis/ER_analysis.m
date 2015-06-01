function allData = ER_analysis()

folders = {'~/data/cohcon/s0300_pilot/s030020150509/'};
ROIs = {'l_v1','l_hmt'};
curAnalysis = {'left_ER'};

tSeries = {};
stimVol = {};
concatInfo = {};

for fi = 1:length(folders)
    [tSeries{end+1},stimVol{end+1},concatInfo{end+1},stimNames] = loadERAnalysis(folder{fi},ROIs,curAnalysis);
end

% At this point each cell is for each folder, and each sub-cell in tSeries
% is for the ROIs. So we can run the analysis by concatenating the
% appropriate runs, and then pulling out of the de-convolution the relevant
% information.

allData = performAnalysis(tSeries,stimVol,concatInfo,stimNames);

function allData = performAnalysis(tSeries,stimVol,concatInfo,stimNames)

for z = 1:length(stimVol) % across folders
    % do concatRuns
end

data = {};

allData = {};
for ri = 1:length(tSeries)
    data{ri} = fitTimecourse(tSeries{ri},stimVol,.5,'fitType=deconv','amplitudeType=fit2');    
    
    % setup info
    [val, cuedTask, valueType] = parseNames(stimNames);
    
    sides = {'left','right'};
    ROIs = {'l_v1','l_hmt'};
    
    [side, roi] = parseROI(ROIs{ri});
    
    %          cued  uncued
    Ns.contrast.cued = [];
    Ns.contrast.uncued = [];
    Ns.coherence.cued = [];
    Ns.coherence.uncued = [];
    contrast.cued = [];
    contrast.uncued = [];
    contrastResp.cued = [];
    contrastResp.uncued = [];
    coherence.cued = [];
    coherenceResp.cued = [];
    coherence.uncued = [];
    coherenceResp.uncued = [];
    
    cueds = {'uncued','cued'};
    values = {'coherence','contrast'};
    resps = {'coherenceResp','contrastResp'};
    for z = 1:length(data{ri}.amplitude)
        evalc(sprintf('%s.%s(end+1) = %f',values{valueType(z)},cueds{(cuedTask(z)==valueType(z))+1},val(z)));
        evalc(sprintf('%s.%s(end+1) = %f',resps{valueType(z)},cueds{(cuedTask(z)==valueType(z))+1},data{ri}.amplitude(z)));
        Ns.(values{valueType(z)}).(cueds{(cuedTask(z)==valueType(z))+1})(end+1) = length(stimVol{z});
    end
    
    for vi = 1:length(values)
        for ci = 1:length(cueds)
            evalc(sprintf('[c,i] = sort(%s.%s)',values{vi},cueds{ci}));
            evalc(sprintf('r = %s.%s(i)',resps{vi},cueds{ci}));
            evalc('curN = Ns.(values{vi}).(cueds{ci})(i)');
            evalc('allData.(sides{side}).(ROIs{roi}).(values{vi}).(cueds{ci}).c = c');
            evalc('allData.(sides{side}).(ROIs{roi}).(values{vi}).(cueds{ci}).r = r');
            evalc('allData.(sides{side}).(ROIs{roi}).(values{vi}).(cueds{ci}).n = curN');
        end
    end
end

function [tSeries, stimVol, concatInfo,stimNames] = loadERAnalysis(folder,allROI,curA)
%%
cdir = pwd;
cd(folder);

%%
view = newView();
view = viewSet(view,'curGroup','Concatenation');
scans = viewGet(view,'nScans');
view = viewSet(view,'curScan',scans);
view = loadAnalysis(view,sprintf('erAnal/%s',curA));
analysis = viewGet(view,'analysis');
d = analysis.d{1};
d.scanNum = 1;
d.groupNum = view.curGroup;
stimNames = d.stimNames;
d = loadroi(d,allROI);
concatInfo = viewGet(view,'concatInfo');
stimVol = d.stimvol;

scanDims = viewGet(view,'scanDims');
r2 = viewGet(view,'overlayData',scans);

tSeries = {};
for ri = 1:length(allROI)
    r = d.roi{roi};
    r.linearScanCoords = sub2ind(scanDims,r.scanCoords(1,:),r.scanCoords(2,:),r.scanCoords(3,:));

    r.r2 = r2(r.linearScanCoords);
    tSeries{end+1} = mean(r.tSeries(r.r2>thresh,:));
end

clear view
clear analysis
clear d

%%

cd(cdir);

function [tSeries] = pullROITSeries(d,roi,thresh)

r = d.roi{roi};
r.linearScanCoords = sub2ind(scanDims,r.scanCoords(1,:),r.scanCoords(2,:),r.scanCoords(3,:));

r.r2 = r2(r.linearScanCoords);
tSeries = mean(r.tSeries(r.r2>thresh,:));


function allData = eventRelatedPedestalPlot(aName)


x = .25:.5:14.75;

    
allData = {};

%% Get the ER analysis
view = newView();
view = viewSet(view,'curGroup','Concatenation');
scans = viewGet(view,'nScans');
view = viewSet(view,'curScan',scans);

view = loadAnalysis(view,'erAnal/left_ER');

analysis = viewGet(view,'analysis');
d = analysis.d{1};
d.r2 = analysis.overlays(1).data{1};



%% Calc Mean Response Across V1? and MT?

d.scanNum = 1;
d.groupNum = view.curGroup;
d = loadroi(d,{'l_v1', 'l_mt'});

%% Scan info
scanDims = viewGet(view,'scanDims');
r2 = viewGet(view,'overlayData',scans);
scm = makescm(view,d.hdrlen,1,d.stimvol);

%% Across ROIS
for roinum = 1:2
    r = d.roi{roinum};
    r.linearScanCoords = sub2ind(scanDims,r.scanCoords(1,:),r.scanCoords(2,:),r.scanCoords(3,:));
    
    r.r2 = r2(r.linearScanCoords);
    
    tSeries = mean(r.tSeries(r.r2>.2,:));
    
    roid = getr2timecourse(tSeries,d.nhdr,d.hdrlen,scm,d.tr);
    
%     figure, hold on
%     for z = 10:2:18
%         errorbar(roid.time,roid.ehdr(z,:),roid.ehdrste(z,:),'color',rand(1,3));
%     end
%     legend(d.stimNames(10:2:18));

    [val, cuedTask, valueType] = parseNames(stimNames);
    [side, roi] = parseROI(r.name);
    
    sides = {'left','right'};
    ROIs = {'v1','mt'};
    
    %          cued  uncued
    Ns.contrast.cued = [];
    Ns.contrast.uncued = [];
    Ns.coherence.cued = [];
    Ns.coherence.uncued = [];
    contrast.cued = [];
    contrast.uncued = [];
    contrastResp.cued = [];
    contrastResp.uncued = [];
    coherence.cued = [];
    coherenceResp.cued = [];
    coherence.uncued = [];
    coherenceResp.uncued = [];
    
    cueds = {'uncued','cued'};
    values = {'coherence','contrast'};
    resps = {'coherenceResp','contrastResp'};
    for z = 1:length(d.amplitude)
        evalc(sprintf('%s.%s(end+1) = %f',values{valueType(z)},cueds{(cuedTask(z)==valueType(z))+1},val(z)));
        evalc(sprintf('%s.%s(end+1) = %f',resps{valueType(z)},cueds{(cuedTask(z)==valueType(z))+1},d.amplitude(z)));
        Ns.(values{valueType(z)}).(cueds{(cuedTask(z)==valueType(z))+1})(end+1) = length(stimVol{z});
    end
    
    for vi = 1:length(values)
        for ci = 1:length(cueds)
            evalc(sprintf('[c,i] = sort(%s.%s)',values{vi},cueds{ci}));
            evalc(sprintf('r = %s.%s(i)',resps{vi},cueds{ci}));
            evalc('curN = Ns.(values{vi}).(cueds{ci})(i)');
            evalc('allData.(sides{side}).(ROIs{roi}).(values{vi}).(cueds{ci}).c = c');
            evalc('allData.(sides{side}).(ROIs{roi}).(values{vi}).(cueds{ci}).r = r');
            evalc('allData.(sides{side}).(ROIs{roi}).(values{vi}).(cueds{ci}).n = curN');
        end
    end
end

clear view

%% Output

% anFolder = '~/proj/att_awe/analysis/csv/';
% if ~exist(anFolder)
%     mkdir(anFolder);
% end
% fname = fullfile(anFolder,'scanData.csv');
% csvwriteh(fname,data,header);

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
%%
val = []; cuedTask = []; valueType = [];

vstr = 'Co';

for i = 1:length(stimNames)
    name = stimNames{i};
    t = str2num(name(strfind(name,'task=')+5:end));
    v = str2num(name(strfind(name,vstr)+4:strfind(name,' and')));
    val(end+1) = v;
    cuedTask(end+1) = t;
    if strfind(name,'Con')
        valueType(end+1) = 2;
    else
        valueType(end+1) = 1;
    end
end


function r = nakarushton(c,n,c50,k,cn)

if ieNotDefined('cn'),cn = c;end

r  = k*(c.^n./(cn.^n+c50.^n));