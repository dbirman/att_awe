function neural = cc_decoding(neural,sid)

%% Step 1: Load Data
folders = neural.folders;
folders = folders(1);

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
    
    ds{fi} = d;
    %
    
    
end

%% Matching Voxels Across Scans
% because the cross-validation is only so good it would be nice to have the
% same voxels across scans--they will be slightly off but it isn't a big
% deal. 

%% Step 2: Get Instances
% now we want to get the average response to every single condition that
% exists in our experiment. We use 'get instances' to do this. Then we want
% to combine this into a giant long-form dataset.
%
adata = [];
aInsts = {};
for fi = 1:length(folders)
    folderz = sprintf('f%s',folders{fi});
    
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
    d = ds{fi};
    
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
    lInst = getInstances(view,lROIs,neural.SCM_f.(folderz).lStim.stimVol,'type','glm','r2cutoff=0.07','hdrlen=15');
    rInst = getInstances(view,rROIs,neural.SCM_f.(folderz).rStim.stimVol,'type','glm','r2cutoff=0.07','hdrlen=15');
    aInst = [lInst rInst];
    
    %%
    tCount = 0;
    for li = 1:length(lInst)
        inst = lInst{li}.classify.instances;
        for ii = 1:length(inst)
            tCount = tCount + length(inst{ii});
        end
    end
    tCount = 0;
    for ri = 1:length(rInst)
        inst = rInst{ri}.classify.instances;
        for ii = 1:length(inst)
            tCount = tCount + length(inst{ii});
        end
    end
    tCount = tCount*100; % max possible size
    
    % setup data to the max possible size, we will remove the end later
    data = zeros(tCount,9);
    count = 1;
    %% calculate left/right data
    disp('Starting left...');
    disppercent(-1/length(lInst));
    for li = 1:length(lInst)
        inst = lInst{li};
        
        if ~isempty(inst.classify.instances)
            sN = neural.SCM_f.(folderz).lStim.stimNames;
            sV = neural.SCM_f.(folderz).lStim.stimVol;
            tsN = neural.SCM.(folderz).lStim.taskNames;
            tsV = neural.SCM.(folderz).lStim.taskSV;
            iV = inst.classify.instanceVol;
            amps = inst.classify.instances;
            [con,coh,task, sv] = sv2long(sV,sN,tsV,tsN);
            
            roi_pos = find(cellfun(@(x) strcmp(x,inst.name),neural.ROIs));
            croi = neural.ROIs{roi_pos};
            [s, r] = parseROI(croi,neural.shortROIs);
            
            
            for ci=1:length(amps)
                if ~isempty(amps)
                    camps = amps{ci};
                    for i = 1:size(camps,1)
                        camp = camps(i,:);
                        vol = iV{ci}(i);
                        
                        sv_pos = find(sv==vol);
                        if ~isempty(sv_pos)
                            ccon = con(sv_pos);
                            ccoh = coh(sv_pos);
                            ctask = task(sv_pos);
                            
                            for ai = 1:length(camp) % for each voxel
                                
                                a = camp(ai);
                                
                                % roi - voxel - con - coh - task - amplitude - instance
                                % - folder - side
                                dat = [r ai ccon ccoh ctask a i fi s];
                                data(count,:) = dat;
                                count = count+1;
                            end
                        end
                    end
                    
                end
            end
        end
        
        disppercent(li/length(lInst));
    end
    disppercent(inf);
    
    disp('Starting right...');
    disppercent(-1/length(rInst));
    for li = 1:length(rInst)
        inst = rInst{li};
        
        if ~isempty(inst.classify.instances)
            
            sN = neural.SCM_f.(folderz).rStim.stimNames;
            sV = neural.SCM_f.(folderz).rStim.stimVol;
            tsN = neural.SCM.(folderz).rStim.taskNames;
            tsV = neural.SCM.(folderz).rStim.taskSV;
            iV = inst.classify.instanceVol;
            amps = inst.classify.instances;
            [con,coh,task, sv] = sv2long(sV,sN,tsV,tsN);
            
            roi_pos = find(cellfun(@(x) strcmp(x,inst.name),neural.ROIs));
            croi = neural.ROIs{roi_pos};
            [s, r] = parseROI(croi,neural.shortROIs);
            
            
            for ci=1:length(amps)
                if ~isempty(amps)
                    camps = amps{ci};
                    for i = 1:size(camps,1)
                        camp = camps(i,:);
                        vol = iV{ci}(i);
                        
                        sv_pos = find(sv==vol);
                        if ~isempty(sv_pos)
                            ccon = con(sv_pos);
                            ccoh = coh(sv_pos);
                            ctask = task(sv_pos);
                            
                            for ai = 1:length(camp) % for each voxel
                                
                                a = camp(ai);
                                
                                % roi - voxel - con - coh - task - amplitude - instance
                                % - folder - side
                                dat = [r ai ccon ccoh ctask a i fi s];
                                data(count,:) = dat;
                                count = count+1;
                            end
                        end
                    end
                    
                end
            end
        end
        disppercent(li/length(rInst));
    end
    disppercent(inf);
    data = data(1:count,:);
    adata = [adata;data];
end

cd(cdir);

%% Write CSV (then can be used for R analyses)

% roi - voxel - con - coh - task - amplitude - instance
%                             % - folder - side
header = {'roi','voxel','con','coh','task','amplitude','instance','folder','side'};
fname = sprintf('~/proj/att_awe/analysis/data/%s_voxels.csv',sid);
fname2 = sprintf('~/proj/Box Sync/cohcon/%s_voxels.csv',sid);
csvwriteh(fname,data,header);
csvwriteh(fname2,data,header);

return

%% Cross-Validated Analysis (Build Encoding Model, Decode from Test Set)

% We are going to do a few things here. First, we split the dataset into 5
% parts *for each voxel in each session*. Then we take the first four folds
