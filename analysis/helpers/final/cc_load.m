function ds = cc_load( neural, sid )
%CC_LOAD Summary of this function goes here
%   Detailed explanation goes here
folders = neural.folders;

if ~isfield(neural,'tSeries')
    neural.tSeries = struct;
end

preF = fullfile('~/data/cohcon/');
cdir = pwd;
name = 'full';
ds = {};
for fi = 1:length(folders)
    %%
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
    
    
    d = loadroi(d,neural.ROIs,'keepNAN=true');
    
    r2 = viewGet(view,'overlayData',d.scanNum);
    
    d.roi = getSortIndex(view,d.roi,r2);
    
    
    % save file
    fname = fullfile(preF,sprintf('%s_loaded_%s.mat',sid,folder));
    save(fname,'d');
    ds{end+1} = fname;
    
end


end

