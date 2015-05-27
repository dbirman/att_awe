function ER_analysis()

cdir = pwd;
cd ~/data/cohcon/s030020150509/

eventRelatedPedestalPlot();
cd(cdir);

function eventRelatedPedestalPlot()

for ti = 1:length(types)

    x = .25:.5:24.75;
    types = {'coherence', 'contrast'};

    data = [];
    %         HDR est  error   con/coh  task   left/right  V1/MT
    header = {'value', 'std', 'pedVal', 'task', 'side', 'roi', 'type', 'voxel'}; 
    type = types{ti};
    %% Get the ER analysis
    view = newView();
    view = viewSet(view,'curGroup','ConcatAll');
    scans = viewGet(view,'nScans');
    view = viewSet(view,'curScan',scans);

    view = loadAnalysis(view,sprintf('erAnal/left_%s_ER',type));

    analysis = viewGet(view,'analysis');
    d = analysis.d{1};
    d.r2 = analysis.overlays(1).data{1};

    %% Calc Mean Response Across V1? and MT?

    d.scanNum = 1;
    d.groupNum = view.curGroup;
    d = loadroi(d,{'left_V1?', 'left_MT?'});

    %% Across ROIS
    for roinum = 1:2
        r = d.roi{roinum};
        roid = [];
        for z = 1:size(r.scanCoords,2)
            for valtask = 1:18
                roid(valtask,z,:) = d.ehdr(r.scanCoords(1,z),r.scanCoords(2,z),r.scanCoords(3,z),valtask,:);
            end
        end

        p1 = find(x>=3.5,1); p2 = find(x>=5,1);

        % v1dsd = squeeze(std(roid,:,p1:p2)); % std across voxels... should be across stimulus presentations, fix later
        roid = squeeze(mean(roid(:,:,p1:p2),3));
    %     roid = squeeze(mean(roid(:,:,p1:p2),1));
    %     roid = squeeze(mean(roid,2));

        [val, task] = parseNames(d.stimNames,type);


        [side, roi] = parseROI(r.name);

        deconvType = strcmp(type,'contrast') + 1;
        for i = 1:size(roid,1)
            for j = 1:size(roid,2)
                data(end+1,:) = [roid(i,j) 0 val(i) task(i) side roi deconvType j];
            end
        end
    end
    
    clear view
end

%% Output

anFolder = '~/proj/att_awe/analysis/csv/';
if ~exist(anFolder)
    mkdir(anFolder);
end
fname = fullfile(anFolder,'scanData.csv');
csvwriteh(fname,data,header);

function [side, roi] = parseROI(roiname)

sname = roiname(1:strfind(roiname,'_')-1);
rname = roiname(strfind(roiname,'_')+1:end-1);

if strcmp(sname,'left')
    side = 1;
else
    side = 2;
end

if strcmp(rname,'V1')
    roi = 1;
else
    roi = 2;
end


function [val, task] = parseNames(names,type)
%%
val = []; task = [];

if strcmp(type,'contrast')
    vstr = 'Con';
else
    vstr = 'Coh';
end

for i = 1:length(names)
    name = names{i};
    t = str2num(name(strfind(name,'task=')+5:end));
    v = str2num(name(strfind(name,vstr)+4:strfind(name,' and')));
    val(end+1) = v;
    task(end+1) = t;
end