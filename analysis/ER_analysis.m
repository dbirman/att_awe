function allData = ER_analysis()

cdir = pwd;
cd ~/data/cohcon/s0300_pilot/s030020150509/

allData = eventRelatedPedestalPlot();
cd(cdir);

function [tSeries, stimVol, concatInfo] = loadERAnalysis(name,roi)
%%
cdir = pwd;
cd(name);

%%
view = newView();
view = viewSet(view,'curGroup','Concatenation');
scans = viewGet(view,'nScans');
view = viewSet(view,'curScan',scans);
view = loadAnalysis(view,'erAnal/left_ER');
analysis = viewGet(view,'analysis');
d = analysis.d{1};
d = loadroi(d,roi);
d.r2 = analysis.overlays(1).data{1};

%%

cd(cdir);

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

    [val, cuedTask, valueType] = parseNames(d.stimNames);
    [side, roi] = parseROI(r.name);
    
    bfit = {};
    for z = 1:size(roid.ehdr,1)
        bfit{z} = fitgamma(roid.ehdr(z,:),x);
        N(z) = length(d.stimvol{z});
    end
    
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
    for z = 1:length(bfit)
        evalc(sprintf('%s.%s(end+1) = %f',values{valueType(z)},cueds{(cuedTask(z)==valueType(z))+1},val(z)));
        evalc(sprintf('%s.%s(end+1) = %f',resps{valueType(z)},cueds{(cuedTask(z)==valueType(z))+1},bfit{z}.params(1)));
        Ns.(values{valueType(z)}).(cueds{(cuedTask(z)==valueType(z))+1})(end+1) = N(z);
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


function [val, cuedTask, valueType] = parseNames(names)
%%
val = []; cuedTask = []; valueType = [];

vstr = 'Co';

for i = 1:length(names)
    name = names{i};
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