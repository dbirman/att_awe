function insts = cc_inst( neural, sid, ds )
preF = fullfile('~/data/cohcon/');

cdir = pwd;
% adata = [];
insts = {};
folders = neural.folders;
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
    
    load(ds{fi});
    
    lindxs = []; rindxs = [];
    for ri = 1:length(d.roi)
        if strcmp(d.roi{ri}.name(1),'l')
            lindxs = [lindxs ri];
        else
            rindxs = [rindxs ri];
        end
    end
    
    lROIs = d.roi(lindxs);
    rROIs = d.roi(rindxs);
    %
    %     lInst = getInstances(view,lROIs,neural.SCM.(folderz).lStim.stimVol,'startLag=6','blockLen=10','minResponseLen=4');
    %     rInst = getInstances(view,rROIs,neural.SCM.(folderz).rStim.stimVol,'startLag=6','blockLen=10','minResponseLen=4');
    lInst = getInstances(view,lROIs,neural.SCM_f.(folderz).rStim.stimVol,'type','glm','r2cutoff=0.05','hdrlen=15','n=inf');
    rInst = getInstances(view,rROIs,neural.SCM_f.(folderz).lStim.stimVol,'type','glm','r2cutoff=0.05','hdrlen=15','n=inf');
    
    fname = fullfile('~/data/cohcon/',sprintf('%s_inst_%s.mat',sid,folder));
    save(fname,'lInst','rInst');
    insts{end+1} = fname;
    
    clear d
end

cd(cdir);
