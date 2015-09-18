function neural = cc_concatER(neural,name)

% Start by concatenating, within each ROI, within each analysis

% Note that the left/right hemisphere stimVols actually correspond to
% opposite mappings (i.e. lCoh or rCoh) of stimuli. But we can combine them
% by ignoring the prefix later on, so it's safe to concatenate everything
% here as long as it's coming form the same ROI.
%%
folders = neural.folders;
shortList = neural.shortROIs;

rois = neural.ROIs;
    
if isfield(neural.tSeries,name)
    neural.tSeries = rmfield(neural.tSeries,name);
end

for ri = 1:length(rois)
            
    roi = rois{ri};
    % get the current ROI that we are in
    [side, rNum] = parseROI(roi,shortList);
    % get the ROI that this maps onto (we just ignore hemisphere)
    rroi = shortList{rNum};
    disp(sprintf('Concatenating for ROI %s to ROI %s tSeries',roi,rroi));
    % Now we do the actual concatentations
    for fi = 1:length(folders)
        folder = folders{fi};
        disp(sprintf('Folder %s data',folder));
        folderz = sprintf('f%s',folder);
        data = neural.tSeries.(folderz).(name);
        sv = neural.SCM.(folderz);

        if isfield(neural.tSeries,name) && isfield(neural.tSeries.(name),rroi)
            % If the concat field already exists, we are going to add to it
            roiTSeries = neural.tSeries.(name).(rroi).tSeries;
            roistimVols = neural.SCM.(name).(rroi).stimVol;
            roiConcat = neural.tSeries.(name).(rroi).concatInfo;
        else
            % If not, we're just going to start from scratch
            roiTSeries = [];
            roistimVols = [];
            roiConcat = {};
        end
        
        cts = data.tSeries{ri};
        
        if any(isnan(cts))
            disp(sprintf('(concatER) Folder %s in ROI %s contains NaNs, not concatenating...',folder,roi));
        else

            if isempty(roiTSeries)
                % If this is the first timeseries we just save the data
                roiTSeries = cts;
                if side==1
                    roistimVols = sv.lStim.(name).stimVol;
                    neural.SCM.(name).(rroi).stimNames = sv.lStim.(name).stimNames;
                else
                    roistimVols = sv.rStim.(name).stimVol;
                    neural.SCM.(name).(rroi).stimNames = sv.rStim.(name).stimNames;
                end
                roiConcat = data.concatInfo;
                % note stimnames is the same across all folders for an ROI
                % so this is okay. BUT: this only applies once prefixes
                % have been ignored obviously...
            else
                % Otherwise we concatenate to the previous data, we also
                % concatenate the stimVols even though the prefixes may be
                % different.
                ccI = data.concatInfo;
                if side==1
                    csv = sv.lStim.(name).stimVol;
                else
                    csv = sv.rStim.(name).stimVol;
                end
                if length(cts)>100
                    [roiTSeries, roistimVols, roiConcat] = concatRuns({roiTSeries, cts}, {roistimVols, csv}, {roiConcat, ccI});
                end
            end
        end
        % Save the data for the next run
        neural.tSeries.(name).(rroi).tSeries = roiTSeries;
        neural.SCM.(name).(rroi).stimVol = roistimVols;
        neural.tSeries.(name).(rroi).concatInfo = roiConcat;
    end
end


function [side, rNum] = parseROI(roi,shortList)

sname = roi(1:strfind(roi,'_')-1);
rname = roi(strfind(roi,'_')+1:end);

if strcmp(sname,'l')
    side = 1;
else
    side = 2;
end

for i = 1:length(shortList)
    if strfind(rname,shortList{i})
        rNum = i;
    end
end