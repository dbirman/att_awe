function SCM = cc_genSCM(neural,sid,group)

%% Navigate Through folders
folders = neural.folders;
SCM = struct;

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
    %     view = loadAnalysis(view,sprintf('erAnal/%s','both_ER'));
    %     analysis = viewGet(view,'analysis');
    %     d = analysis.d{neural.concatScan};
    %     d.scanNum = neural.concatScan;
    %     d.groupNum = view.curGroup;
    
    if isempty(group)
        prefixes = {'r','l'};
        scmGroups = {'lStim','rStim'};
        for si = 1:length(scmGroups)
            scmGroup = scmGroups{si};
            prefix = prefixes{si};
            
            allStims = {sprintf('%sCon_x_%sCoh',prefix,prefix)};
            
            if ~isfield(SCM,folderz)
                SCM.(folderz) = struct;
            end
            if ~isfield(SCM.(folderz),scmGroup)
                SCM.(folderz).(scmGroup) = struct;
            end
            [SCM.(folderz).(scmGroup).stimVol, SCM.(folderz).(scmGroup).stimNames, ~] = getStimvol(view,allStims);
            [SCM.(folderz).(scmGroup).taskSV, SCM.(folderz).(scmGroup).taskNames, ~] = getStimvol(view,'nTask');
        end
    else
        [SCM.(folderz).stimVol, SCM.(folderz).stimNames, ~] = getStimvol(view,group);
    end
end
cd(cdir);