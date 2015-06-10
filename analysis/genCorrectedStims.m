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

%% Get current stimulus values

randVars = task{1}{1}.randVars;
trials = task{1}{1}.trialnum;
param = task{1}{1}.block.parameter;

if ~isfield(randVars,'coherence')
    disp('This stimfile doesn''t have a coherence field, skipping.');
    return
end

baseCoh = randVars.coherence(1:trials);
baseCon = randVars.contrast(1:trials);

lCon = []; lCoh = [];
rCon = []; rCoh = [];

if any(param.catch(1:trials)>0)
    nTask = randVars.task(1:trials) + param.catch(1:trials) .* 2;
else
    nTask = randVars.task(1:trials);
end

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
lCon = nearestValue(lCon,2); rCon = nearestValue(rCon,2);
lCoh = nearestValue(lCoh,1); rCoh = nearestValue(rCoh,1);

if ~isfield(randVars,'lCon'), addCalculatedVar('lCon',lCon,file,'backup=0'); end
if ~isfield(randVars,'lCoh'), addCalculatedVar('lCoh',lCoh,file,'backup=0'); end
if ~isfield(randVars,'rCon'), addCalculatedVar('rCon',rCon,file,'backup=0'); end
if ~isfield(randVars,'rCoh'), addCalculatedVar('rCoh',rCoh,file,'backup=0'); end
if ~isfield(randVars,'nTask'), addCalculatedVar('nTask',nTask,file,'backup=0'); end

function val = nearestValue(val,type)

if length(val)>1
    for i = 1:length(val)
        val(i) = nearestValue(val(i),type);
    end
    return
end

if type==1
    % coherence
    vList = [0 .1 .25 .7];
else
    % contrast
    vList = [.2 .4 .6 .8];
end

if val > 1 || val < 0
    disp('only works on percentages');
    keyboard
end

high = find(vList >= val,1);
low = find(vList <= val,1);

if abs(vList(high)-val) < abs(vList(low)-val)
    val = vList(high);
else
    val = vList(low);
end