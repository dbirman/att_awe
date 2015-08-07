function genCorrectedStims(folder)

%% If isempty, default to 

%% Get stim files from folder

files = dir(fullfile(folder,'Etc','*.mat'));

%% Mkdir

if ~isdir(fullfile(folder,'Etc','Orig'))
    mkdir(fullfile(folder,'Etc','Orig'));
end

%% Run for each
for i = 1:length(files)
    modifyStim(fullfile(folder,'Etc',files(i).name),fullfile(folder,'Etc','Orig',files(i).name));
end

function modifyStim(file,backup)

%% Load stim

load(file);

copyfile(file,backup);

%% Track file loc

%% Check for vals

if ~exist('stimulus','var') || ~exist('task','var') || ~exist('myscreen','var')
    disp('Load failed for this stimfile, check by hand?');
    keyboard
end

%% nTask

randVars = task{1}{1}.randVars;
if ~isfield(randVars,'coherence')
    disp('This stimfile doesn''t have a coherence field, skipping.');
    return
end



trials = task{1}{1}.trialnum;
param = task{1}{1}.block.parameter;
if any(param.catch(1:trials)>0)
    nTask = randVars.task(1:trials) + param.catch(1:trials) .* 2;
else
    nTask = randVars.task(1:trials);
end

if ~isfield(randVars,'nTask')
    addCalculatedVar('nTask',nTask,file,'backup=0');
end

%% xCon and xCoh


baseCoh = randVars.coherence(1:trials);
baseCon = randVars.contrast(1:trials);

lCon = []; lCoh = [];
rCon = []; rCoh = [];
if ~isfield(randVars,'lCon')
    for t = 1:trials
        if param.conSide(t) == 1
            % contrast left
            lCon(t) = baseCon(t) + randVars.conDelta(t);
            rCon(t) = baseCon(t);
        else
            lCon(t) = baseCon(t);
            rCon(t) = baseCon(t) + randVars.conDelta(t);
        end
        if param.cohSide(t) == 1
            lCoh(t) = baseCoh(t) + randVars.cohDelta(t);
            rCoh(t) = baseCoh(t);
        else
            lCoh(t) = baseCoh(t);
            rCoh(t) = baseCoh(t) + randVars.cohDelta(t);
        end
    end
    addCalculatedVar('lCon',lCon,file,'backup=0');
    addCalculatedVar('rCon',rCon,file,'backup=0');
    addCalculatedVar('lCoh',lCoh,file,'backup=0');
    addCalculatedVar('rCoh',rCoh,file,'backup=0');
end
