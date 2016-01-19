function dsfiles = cc_buildset( neural, sid, insts, rois )

dsfiles = {};

folders = neural.folders;
for fi = 1:length(folders)
    
    folder = folders{fi};
    folderz = sprintf('f%s',folder);
    
    load(insts{fi}); % this loads lInst and rInst
    
    %%
    % restrict to specific ROIs
    r_idxs = [];
    n_voxsl = 0; n_voxsr = 0;
    for ri = 1:length(rois)
        roi = rois{ri};
        for rii = 1:length(neural.shortROIs)
            if strcmp(roi,neural.shortROIs{rii})
                r_idxs = [r_idxs rii];
                n_voxsl = n_voxsl + lInst{rii}.n;
                n_voxsr = n_voxsr + rInst{rii}.n;
            end
        end
    end
    
    %%l
    lInst = lInst(r_idxs); % drop all rois we don't need
    rInst = rInst(r_idxs);
    
    %% Now, we build our dataset:
    
    % format:
    % con - coh - task - side - vox1 vox2 ... voxN

    %%
    ldata = zeros(100000,4+n_voxsl);
    rdata = zeros(100000,4+n_voxsr);
    count = 1;
    
    %% calculate left/right data
     
    % pull out stimvols, we will use to parse the data
    all_sv = lInst{1}.classify.stimvol;
%     all_iv = lInst{1}.classify.instanceVol; % should be identical to sv?
    % pull out the task data
    sN = neural.SCM_f.(folderz).rStim.stimNames; % right stim goes with left voxels
    sV = neural.SCM_f.(folderz).rStim.stimVol;
    tsN = neural.SCM.(folderz).rStim.taskNames;
    tsV = neural.SCM.(folderz).rStim.taskSV;
    [con,coh,task, sv] = sv2long(sV,sN,tsV,tsN);
    
    %%
    % for each stimvol set
    disppercent(-1/length(all_sv),'LEFT');
    for si = 1:length(all_sv)
        c_sv = all_sv{si};
        
        % for each stimvol
        for ci = 1:length(c_sv)
            cur = c_sv(ci);
            
            pos = find(cur==sv);
            ccon = con(pos); ccoh = coh(pos); ctask = task(pos);
            
            if ~isempty(pos)
                dat = [ccon ccoh ctask 1];
                % for each roi
                for ri = 1:length(r_idxs)
                    amps = lInst{ri}.classify.instances{si}(ci,:);
                    dat = [dat amps];
                end
                ldata(count,:) = dat;
                count = count + 1;
            end
        end    
        disppercent((si-1)/length(all_sv));
    end
    disppercent(inf);
    ldata = ldata(1:count,:);
    
    %% calculate left/right data
     
    % pull out stimvols, we will use to parse the data
    all_sv = rInst{1}.classify.stimvol;
%     all_iv = lInst{1}.classify.instanceVol; % should be identical to sv?
    % pull out the task data
    sN = neural.SCM_f.(folderz).lStim.stimNames; % right stim goes with left voxels
    sV = neural.SCM_f.(folderz).lStim.stimVol;
    tsN = neural.SCM.(folderz).lStim.taskNames;
    tsV = neural.SCM.(folderz).lStim.taskSV;
    [con,coh,task, sv] = sv2long(sV,sN,tsV,tsN);
        %%
    % for each stimvol set
    disppercent(-1/length(all_sv),'RIGHT');
    for si = 1:length(all_sv)
        c_sv = all_sv{si};
        
        % for each stimvol
        for ci = 1:length(c_sv)
            cur = c_sv(ci);
            
            pos = find(cur==sv);
            ccon = con(pos); ccoh = coh(pos); ctask = task(pos);
            
            if ~isempty(pos)
                dat = [ccon ccoh ctask 2];
                % for each roi
                for ri = 1:length(r_idxs)
                    amps = rInst{ri}.classify.instances{si}(ci,:);
                    dat = [dat amps];
                end
                rdata(count,:) = dat;
                count = count + 1;
            end
        end    
        disppercent((si-1)/length(all_sv));
    end
    disppercent(inf);
    rdata = rdata(1:count,:);
    
    fname = fullfile('~/data/cohcon/',sprintf('%s_dataset_%s.mat',sid,folder));
    dsfiles{end+1} = fname;
    save(fname,'ldata','rdata');
end

