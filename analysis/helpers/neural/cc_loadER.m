function neural = cc_loadER(neural,sid,name)

%% Navigate Through folders
folders = neural.folders;

if ~isfield(neural,'tSeries')
    neural.tSeries = struct;
end

preF = fullfile('~/data/cohcon/');
cdir = pwd;
for fi = 1:length(folders)
    folder = folders{fi};
    folderz = sprintf('f%s',folder);
    fullFolder = fullfile(preF,sprintf('%s%s',sid,folder));
    
    mrQuit;
    cd(fullFolder);
    
    %% Setup a view, don't load ROIs
    view = newView();
    view = viewSet(view,'curGroup','Concatenation');
    view = viewSet(view,'curScan',neural.concatScan);
    view = loadAnalysis(view,sprintf('erAnal/%s','both_ER'));
    analysis = viewGet(view,'analysis');
    d = analysis.d{neural.concatScan};
    d.scanNum = neural.concatScan;
    d.groupNum = view.curGroup;
    thresh = neural.r2_cutoffs(fi);
    d = loadroi(d,neural.ROIs);
%     d = erBootstrap('view',view,'d',d,'iterations=10');

    concatInfo = viewGet(view,'concatInfo');

    tSeries = {}; % meaned tSeries
%     fTSeries = {}; % all voxels surviving r^2 cutoff

    scanDims = viewGet(view,'scanDims');
    r2 = viewGet(view,'overlayData',d.scanNum);
    for ri = 1:length(d.roi)
        r = d.roi{ri};
        r.linearScanCoords = sub2ind(scanDims,r.scanCoords(1,:),r.scanCoords(2,:),r.scanCoords(3,:));

        r.r2 = r2(r.linearScanCoords);
%         fTSeries{end+1} = r.tSeries(r.r2>thresh,:);
        tSeries{end+1} = nanmean(r.tSeries(r.r2>thresh,:));
    end

%     neural.tSeries.(folderz).(name).fullTSeries = fTSeries;
    neural.tSeries.(folderz).(name).ROIs = neural.ROIs;
    neural.tSeries.(folderz).(name).tSeries = tSeries;
    neural.tSeries.(folderz).(name).concatInfo = concatInfo;
end

%%
cd(cdir);
